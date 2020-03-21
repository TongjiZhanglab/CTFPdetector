
####
####  0. Configure package
####
'install package'
if(!require(data.table)) {
    install.packages('data.table', dep=T)
}
if(!require(ggplot2)) {
    install.packages('ggplot2', dep=T)
}
if(!require(RColorBrewer)) {
    install.packages('RColorBrewer', dep=T)
}
if(!require(testthat)) {
    install.packages('testthat', dep=T)
}
if(!require(FME)) {
    install.packages('FME', dep=T)
}
if(!require(parallel)) {
    install.packages('parallel', dep=T)
}
if(!require(DESeq2)){
    if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
    }
    BiocManager::install("DESeq2")
}

'Load package'
library(testthat)
library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(FME)
library(parallel)


####
####  0. Data preprocessing
####

'0. tools'
rm.dup <- function(path, sep){
    rawdata <- read.csv(path, header=T, sep=sep, as.is=T)
    X <- as.data.table(rawdata)
    keys = colnames(rawdata)[1]
    data <- as.data.frame(X[,lapply(.SD,mean),keys])
    rownames(data) <- data[,1]
    data <- data[,-1]
    return(data)
}


'1. Data split(from dataframe to list)'
split.data.bytime <- function(data){
    print('Split data by time point...')
    data.list <- list()
    time.name <- unique(as.character(data[1, ]))
    time.cumsum <- cumsum(table(as.character(data[1, ]))[time.name])
    data <- data[-1, ]
    for (i in 1:length(time.cumsum)){
        if (i == 1){
            data.list[[time.name[i]]] <- data.matrix(data[,seq(time.cumsum[i])])
        }else{
            data.list[[time.name[i]]] <- data.matrix(data[,seq(time.cumsum[i-1]+1, time.cumsum[i])])
        }
    }
    return(data.list)
}


'2. Data quality control'
savefig.ggplot <- function(fig, name, width=8, height = 6){
    pdf(name, width=width, height=height)
    print(fig)
    dev.off()
}

filter.cells.by.MT.pect <- function(data, path, MTcutoff=0.2){
    print('Filter cells by mitochondria reads percentage...')
    mt.pect <- c() # 存储对应每个细胞的线粒体reads占比
    MTidx <- grep("^MT-", rownames(data[[1]]))
    for (i in names(data)){
        MTexp <- colSums(data[[i]][MTidx, ])
        dataexp <- colSums(data[[i]])
        MT.pect <- MTexp/dataexp
        mt.pect <- rbind(mt.pect, data.frame(time=rep(i, length(dataexp)), MT.pect=MT.pect*100))
        keepcol <- MT.pect<MTcutoff
        datakeep <- data[[i]][,keepcol]
        data[[i]] <- datakeep
        print(sprintf("Filter cells by mitochondria genes' reads percentage. Keep %s%% cells (total cells: %s)", 
          as.character(formatC(sum(keepcol)/length(keepcol)*100), 2), as.character(length(keepcol))))
    }
    print('Quality control plot...')
    p<-ggplot(mt.pect, aes(x=time, y=MT.pect, fill=time)) +
      geom_violin(trim=FALSE) + labs(title="Quality Control ") + 
      geom_hline(yintercept = c(20))
    savefig.ggplot(p, paste(path, 'QC.pdf', sep='/'), width=8, height=6)
    return(data)
}


'3. Data normalization'
get.sizefactor <- function(data){
    coldata <- data.frame(condition = colnames(data))
    rownames(coldata) <- colnames(data)
    dds <- DESeqDataSetFromMatrix(countData = round(data), colData = coldata, design = ~ condition)
    dds <- estimateSizeFactors(dds)
    size.factor <- sizeFactors(dds)
    return(size.factor)
}

get.norm.data <- function(data, path){
    print('Data normalization...')
    read.cnt <- c()
    for (i in names(data)){
        count.depth <- colSums(data[[i]])
        size.factor <- get.sizefactor(data[[i]])
        data[[i]] <- sweep(data[[i]], 2, size.factor, `/`) # normalize by size factor
        min.val <- min(colSums(data[[i]]))
        data[[i]] <- data[[i]]/min.val*1000000 # transform to million scale
        data[[i]] <- log(data[[i]]+1, 2) # log2 transform
        rownames(data[[i]]) <- toupper(rownames(data[[i]]))
        read.cnt <- rbind(read.cnt, data.frame(count_depth=count.depth, size_factor=size.factor))
    }
    pdf(paste(path, 'Normalization.pdf', sep='/'), width=8, height=6)
    plot(read.cnt[,'size_factor'] , read.cnt[,'count_depth'],  main="Normalization" , pch=20 , cex=1 , col=rgb(0.3,0.5,1,0.4)  , xlab="size.factor" , ylab="count.depth" )
    dev.off()
    return(data)
}


'4. Keep more variable TF'
get.TF.exp <- function(data, TF){
    overlapTF <- intersect(rownames(data[[1]]), TF)
    for (i in names(data)){
        data[[i]] <- data[[i]][overlapTF, ]
    }
    return(data)
}

get.variable.TF <- function(data, TFs, path, min.cutoff=0.1){
	data <- get.TF.exp(data, TFs) # 找出TF表达
    alldata <- c()
    for (i in names(data)){
        alldata <- cbind(alldata, data[[i]])
    }
    gene.mean <- apply(alldata, 1, mean)
    gene.sd <- apply(alldata, 1, sd)
    CV <- gene.sd/(gene.mean+mean(gene.mean))
    var.idx <- CV>min.cutoff
    colors <- c()
    for (i in as.numeric(var.idx)){
        if (i==1){
            colors <- c(colors, 'pink')
        }else{
            colors <- c(colors, 'grey')
        }
    }
    for (i in names(data)){
        data[[i]] <- data[[i]][var.idx, ]
    }
    pdf(paste(path, 'variableTF.pdf', sep='/'), width=8, height=6)
    plot(gene.mean, gene.sd, col=colors, main="Variable TF Filter" , pch=20 , cex=1 , xlab="gene.mean" , ylab="gene.sd")
    dev.off()
    return(data)
}




####
####  1. Core TF filter
####

#(1) 数据离散化
get_discrate_data <- function(TF_Days){
  TF_B=list()
  for(i in 1:length(TF_Days)){
    tf_days=TF_Days[[i]]
    tf_B=c()
    for(k in 1:ncol(tf_days)){
      ex=tf_days[,k]
      ex[ex>0.1]=1
      ex[ex<=0.1]=0
      tf_B=cbind(tf_B,ex)
    }
    colnames(tf_B)=colnames(tf_days)
    rownames(tf_B)=rownames(tf_days)
    TF_B[[i]]=tf_B
  }
  names(TF_B)=names(TF_Days)
  return(TF_B)
}

#(2) 构建每个细胞内部的邻接矩阵
get_Network_Adj <- function(TF_B){
  Network_Adj=list()
  for(i in 1:length(TF_B)){
    tf_B=TF_B[[i]]
    network_adj=list()
    for(cell in 1:ncol(tf_B) ){
      #intial_matrix=matrix(1,nrow = nrow(tf_B),ncol=nrow(tf_B))-diag(rep(1,nrow(tf_B)))
      d=tf_B[,cell]
      #interaction=d%*%t(d)
      #network_adj_cell=interaction*intial_matrix
      network_adj_cell=d%*%t(d)
      rownames(network_adj_cell)=colnames(network_adj_cell)
      network_adj[[cell]]=network_adj_cell
    }
    names(network_adj)=colnames(tf_B)
    Network_Adj[[i]]=network_adj
  }
  
  names(Network_Adj)=names(TF_B)
  return(Network_Adj)
}

#(3)整合state的邻接矩阵
get_state_Network_Adj <- function(Network_Adj){
  
  Sum_adj=list()
  for(i in 1:length(Network_Adj)){
    
    network_adj=Network_Adj[[i]]
    sum_adj=0
    for (k in 1:length(network_adj)){
      adj=network_adj[[k]]
      sum_adj=sum_adj+adj
    }
    
    
    Sum_adj[[i]]=sum_adj
  }
  names(Sum_adj)=names(Network_Adj)
  return(Sum_adj)
}

#(4)定义一个熵
test.entropy <- function(d){
  res <- 0
  for(i in 1:length(d))
  {
    if(d[i]!=0)
      res <- res + d[i]*log2(d[i])
  }
  return (-res)
}

get_Edge_ENTROPY <- function(Integrate_Network_Adj,Ncell){
  
  ENTROPY=c()
  for(k in 1:length(Integrate_Network_Adj)){
    integrate_Network_Adj=Integrate_Network_Adj[[k]]
    Entropy=c()
    ncell=Ncell[[k]]
    for(i in 1:ncol(integrate_Network_Adj)){
      
      gene=colnames(integrate_Network_Adj)[i]
      adj=integrate_Network_Adj[-i,i]
      vec=integrate_Network_Adj[i,i]
      if(vec> 0.05*ncell){
        com_prob=adj/ncell
        prob=com_prob/sum(com_prob)
        
        entropy=(vec/ncell)*test.entropy(prob)
      } else {
        entropy=0
      }
      
      Entropy=append(Entropy,entropy)}
    names(Entropy)=colnames(integrate_Network_Adj)
    ENTROPY=cbind(ENTROPY,Entropy)
  }
  colnames(ENTROPY)=names(Integrate_Network_Adj)
  return(ENTROPY)
}

#(5)判断是否为符合标准的FDE
get_Judgement <- function(ENTROPY){
  
  nstate=ncol(ENTROPY)
  Judgement=c()
  for(i in 1:nrow(ENTROPY)){
    entropy_gene=ENTROPY[i,]
  
    sign_diff=sign(entropy_gene[2:nstate]-entropy_gene[1:(nstate-1)])

    state_diff_label=names(sign_diff)
    judgement="none"
    for(k in 1:(length(state_diff_label)-1)){
      if ((sum(sign_diff[1:k]==1)+sum(sign_diff[1:k]==0))==k   & sum(sign_diff[(k+1):length(state_diff_label)]== -1)==(length(state_diff_label)-k) ){
        judgement = state_diff_label[k]
      }
    }

    cur_cnt = 0
    max_cnt = 0
    for (k in 1:length(sign_diff)){
        if (sign_diff[k]==1){
            cur_cnt = 0
            if (cur_cnt > max_cnt){
                max_cnt = cur_cnt
            }
        }else if (sign_diff[k]==(-1) | sign_diff[k]==0){
            cur_cnt = cur_cnt + 1
            if (cur_cnt > max_cnt){
                max_cnt = cur_cnt
            }

        }
    }
    if (max_cnt>=(length(sign_diff)-1)){
      judgement = state_diff_label[1]
    }
    Judgement=append(Judgement, judgement)
}
  names(Judgement)=rownames(ENTROPY)
  return(Judgement)
}

#(6)对三类核心转录因子标识物进行筛选
CTF_filtering <- function(expRdata){
  # 将数据离散化
  TF_B <- get_discrate_data(expRdata)
  
  # 构建每个细胞内部的邻接矩阵
  Network_Adj <- get_Network_Adj(TF_B)
  
  # 整合state的邻接矩阵
  Integrate_Network_Adj <- get_state_Network_Adj(Network_Adj)
  
  # 定义一个熵
  Ncell <- lapply(Network_Adj,length)
  ENTROPY <- get_Edge_ENTROPY(Integrate_Network_Adj,Ncell)
  
  # 判断FDE
  Judgement <- get_Judgement(ENTROPY)
  keepTFs <- names(Judgement[!Judgement=='none'])

  for (i in names(expRdata)){
    expRdata[[i]] <- expRdata[[i]][keepTFs, ]
  }
  return(expRdata)
}




####
####  2. Core TF rank
####
'1. Find core cell fate TFs'
'TF and expression gene overlap'
# 计算表达list和TF的gene name的overlap
gene_overlap.F <- function(explist,tfpath){
    TF_names <- toupper(as.character(read.csv(tfpath,sep='\t',header=F,as.is=T)[,1]))
    name_list <- list()
    for (i in 1:length(explist)){
        name_list[[i]] <- toupper(as.character(rownames(explist[[i]])))
    }
    name_list[[length(explist)+1]] <- TF_names
    overlap_genes <- Reduce(intersect, name_list)
    return(overlap_genes)
}


'Calculate one gene fitness of sigmoid function'
# 计算某个基因的sigmoid曲线的拟合优度和参数以及稳定状态的时间点
get_one_gene_logistic_parameter.F <- function(exp_list,gene,state_start,state_end){
    # 计算每个时间点的样本数
    num_vec <- c(0)
    num_add <- 0
    for (i in 1:length(exp_list)){
        num <- ncol(exp_list[[i]])
        num_add <- num_add+num
        num_vec <- c(num_vec,num_add)
    }
    # 某基因表达值
    exp <- c()
    for (i in state_start:state_end){
        exp <- c(exp,as.numeric(exp_list[[i]][gene,])) 
        # exp <- c(exp,sort(as.numeric(exp_list[[i]][gene,]))) #################
    }
    exp[exp==0] <- 0.001
    data <- data.frame('x'=seq(length(exp)),'y'=exp)
    # 拟合
    NM=try(nls_model <- nls(y~ SSlogis(x, Asym, xmid, scal), data),silent = T) #Asym渐近线，xmid 1/2渐近线对应的x值，scale对x轴缩放
    # print(NM)
    if(!"try-error" %in% class(NM)){
        # print(gene)
        goodness_of_fit <- cor(as.numeric(data[,'y']),predict(nls_model)) #拟合优度
        coefficient <- coef(nls_model) #sigmoid曲线估计出来的3个参数
        pre <- predict(nls_model) #sigmoid曲线预测值
        if(coefficient[3]>0){ #也就是scale
            vari_region_high <- max(which(coefficient[1]-pre >0.1))
            vari_region_low <- min(which( pre >0.1))
        }else{
            vari_region_high <- min(which(coefficient[1]-pre >0.1))
            vari_region_low <- max(which( pre >0.1))
        }
        for (k in 1:(length(num_vec)-1)){
            if (vari_region_high+num_vec[state_start]>num_vec[k] & vari_region_high+num_vec[state_start]<=num_vec[k+1]){
                vari_region_high <- k
            }
            if (vari_region_low+num_vec[state_start]>num_vec[k] & vari_region_low+num_vec[state_start]<=num_vec[k+1]){
                vari_region_low <- k
            }
        }
        out <- t(data.frame(c(goodness_of_fit,coefficient,names(exp_list)[vari_region_low],names(exp_list)[vari_region_high])))
    }else{
        out <- t(data.frame(rep(NA,6)))
    }
    # plot(as.numeric(data[,'x']),pre,type='l',col='red',main=gene)
    rownames(out) <- gene
    colnames(out) <- c("goodness_of_fit","Asym","xmid","scale","stablesite_low","stablesite_high")
    return(out)
}


'Calculate all genes fitness of sigmoid function(fill NA if can not fit)'
# 计算多个基因的sigmoid拟合优度和参数以及稳定状态,如果某gene无法拟合则填充NA
get_genes_logistic_parameter.F <- function(exp_list,genes,state_start,state_end){
    out_dataframe <- c()
    for (i in genes){
        log_parameter <- get_one_gene_logistic_parameter.F(exp_list,i,state_start,state_end)
        out_dataframe <- rbind(out_dataframe,log_parameter)
    }
    return(out_dataframe)
}


'Calculate all fitness/scale from start time to all other timepoint'
# 计算所有时间点对start时间点过程的fit_goodness和scale
all_time_fitgoodness.scale.F <- function(exp_list,genes){
    all_scale <- c()
    all_fit_goodness <- c()
    for (i in 2:length(exp_list)){
        cur.time.fit <- get_genes_logistic_parameter.F(exp_list,genes,1,i)
        cur.fit_goodness <- as.numeric(cur.time.fit[,1])
        cur.scale <- as.numeric(cur.time.fit[,4])
        all_fit_goodness <- cbind(all_fit_goodness, cur.fit_goodness)
        all_scale <- cbind(all_scale, cur.scale)
    }
    rownames(all_fit_goodness) <- genes;rownames(all_scale) <- genes
    fitness.col <- c();for (i in 2:length(exp_list)){fitness.col <- c(fitness.col,paste('fit_goodness',names(exp_list)[i],sep='_'))}
    scale.col <- c();for (i in 2:length(exp_list)){scale.col <- c(scale.col,paste('scale',names(exp_list)[i],sep='_'))}
    colnames(all_fit_goodness) <- fitness.col;colnames(all_scale) <- scale.col
    return(list(all_fit_goodness, all_scale))
}


'Calculate one gene t-test statistic value between predicted and true expression'
# 计算某个基因的差异统计量ttest
get_gene_t_statistic <- function(exp_list,all_logistic_fitness,all_logistic_scale){
    goodness_of_fits <- c()
    scales <- c()
    stats <- c()
    best_time_points <- c()
    for (i in 1:nrow(all_logistic_fitness)){
        whichmax <- as.numeric(which.max(all_logistic_fitness[i,]))+1
        if (length(whichmax)==0){
            '如果各个时间点对起始点的fitness都是NA,则都填充NA'
            goodness_of_fit <- NA
            scale <- NA
            stat <- NA
            best_time_point <- NA
        }else{
            if (whichmax<length(exp_list)){
                '拟合度最高的时间点到最开始时间点的表达值'
                exp <- c()
                for (j in 1:whichmax){
                    exp <- c(exp,as.numeric(exp_list[[j]][rownames(all_logistic_fitness)[i],]))
                    # exp <- c(exp,sort(as.numeric(exp_list[[j]][rownames(all_logistic_fitness)[i],]))) ################
                }
                exp[exp==0] <- 0.001
                exp_before <- data.frame('x'=seq(length(exp)),'y'=exp)
                '拟合度最高的时间点到最后一个时间点的表达值'
                exp <- c()
                for (k in (whichmax+1):length(exp_list)){
                    exp <- c(exp,as.numeric(exp_list[[k]][rownames(all_logistic_fitness)[i],]))
                }
                exp[exp==0] <- 0.001
                exp_after <- data.frame('x'=seq((nrow(exp_before)+1),(nrow(exp_before)+length(exp))),'y'=exp)
                '用exp before数据进行拟合sigmoid曲线'
                nn <- try(fm1DNase1 <- nls(y~ SSlogis(x, Asym, xmid, scal), exp_before),silent = T)
                # print(rownames(all_logistic_fitness)[i])
                goodness_of_fit <- cor(as.numeric(exp_before$y),predict(fm1DNase1))
                coefficient <- coef(fm1DNase1)
                pre_after <- predict(fm1DNase1,data.frame(x=exp_after[,'x']))
                if (coefficient[3]>0){
                    try_ttest <- try(t <- t.test(pre_after-as.numeric(exp_after[,'y']),alternative = "greater"),silent = T)
                    if(!'try-error' %in% class(try_ttest)){
                        stat=as.numeric(t$statistic)
                    }else{
                        stat=NA
                    }
                }else{
                    try_ttest <- try(t <- t.test(as.numeric(exp_after[,'y'])-pre_after,alternative = "greater"),silent = T)
                    if(!'try-error' %in% class(try_ttest)){
                        stat <- as.numeric(t$statistic)
                    }else{
                        stat <- NA
                    }
                }
                scale <- as.numeric(coefficient[3])
                best_time_point <- names(exp_list)[whichmax]
            }else{
                '如果fit最好的是最后一个时间点则t statistics为NA其余都有值'
                stat <- NA
                goodness_of_fit <- all_logistic_fitness[i,ncol(all_logistic_fitness)]
                scale <- as.numeric(all_logistic_scale[i,ncol(all_logistic_scale)])
                best_time_point <- names(exp_list)[length(exp_list)]
            }
        }
        goodness_of_fits <- c(goodness_of_fits,goodness_of_fit)
        scales <- c(scales,scale)
        stats <- c(stats,stat)
        best_time_points <- c(best_time_points,best_time_point)
    }
    out <- cbind(goodness_of_fits,scales,stats,best_time_points)
    colnames(out) <- c('goodness_of_fit_best','scale','t_statistic','best_time_points')
    rownames(out) <- rownames(all_logistic_fitness)
    out.filter <- out[!is.na(as.numeric(out[,3])),] # 移除最后一个时间点fit最好的基因
    return(out.filter)
}


'core-transcription-factor finding'
CTF_finding <- function(expRdata, TFpath, path){
    'load data'
    for (i in 1:length(expRdata)){
        rownames(expRdata[[i]]) <- toupper(rownames(expRdata[[i]]))
    }
    Days <- expRdata
    if (!is.null(TFpath)){
        TF_genes <- gene_overlap.F(Days,TFpath) #TF genes
        'TF expression list'
        for (i in 1:length(Days)){
            Days[[i]] <- Days[[i]][TF_genes,] #TF expression
        }
    }else{
        allnames <- list()
        for (j in 1:length(Days)){
            allnames[[j]] <- toupper(rownames(Days[[j]]))
        }
        TF_genes <- Reduce(intersect, allnames)
    }

    'up and down gene'
    all <- get_genes_logistic_parameter.F(Days,TF_genes,1,length(Days))
    all <- all[!apply(all, 1, function(x) any(is.na(x))),] # filter NA line
    all_up <- all[as.numeric(all[,4])>0, ] # scale >0
    all_down <- all[as.numeric(all[,4])<0, ] # scale <0

    'all time fitness and scale'
    logistic_fitness_scale <- all_time_fitgoodness.scale.F(Days,TF_genes)
    all_t_statistic <- get_gene_t_statistic(Days, logistic_fitness_scale[[1]], logistic_fitness_scale[[2]])
    all_t_statistic_up <- all_t_statistic[as.numeric(all_t_statistic[,2])>0,] # fitness line is from low to high

    'filter updown, remove the median expression value of best point is 0'
    rm.idx <- c()
    for (i in rownames(all_t_statistic_up)){
        time <- all_t_statistic_up[i, 'best_time_points']
        gene.med <- median(as.numeric(expRdata[[time]][i,]))
        rm.idx <- c(rm.idx, gene.med>0)
    }
    all_t_statistic_up <- all_t_statistic_up[rm.idx, ]

    'Plot fitness and t-statistics distribution'
    pdf(paste(path, 'Fitness.tstatistic.distribution.pdf', sep='/'))
    # print(path)
    par(mfrow=c(2,2))
    plot(sort(as.numeric(all_up[, 1])), seq(nrow(all_up))/nrow(all_up), col='#636363', pch=20, type = "l", lwd=4, 
        xlab='Fitness', ylab=paste('Cumulative Ratio (total gene: ', nrow(all_up), ')', sep=''), main='Distribution of fitness(up)')
    plot(sort(as.numeric(all_down[, 1])), seq(nrow(all_down))/nrow(all_down), col='#636363', pch=20, type = "l", lwd=4, 
        xlab='Fitness', ylab=paste('Cumulative Ratio (total gene: ', nrow(all_down), ')', sep=''), main='Distribution of fitness(down)')
    plot(sort(as.numeric(all_t_statistic_up[, 1])), seq(nrow(all_t_statistic_up))/nrow(all_t_statistic_up), col='#636363', pch=20, type = "l", lwd=4, 
        xlab='Fitness', ylab=paste('Cumulative Ratio (total gene: ', nrow(all_t_statistic_up), ')', sep=''), main='Distribution of fitness(updown)')
    plot(sort(as.numeric(all_t_statistic_up[, 3])), seq(nrow(all_t_statistic_up))/nrow(all_t_statistic_up), col='#636363', pch=20, type = "l", lwd=4, 
        xlab='t-statistics', ylab=paste('Cumulative Ratio (total gene: ', nrow(all_t_statistic_up), ')', sep=''), main='Distribution of t statistics(updown)')
    dev.off()

    'return data'
    out <- list(all_up, all_down, all_t_statistic_up, TF_genes)
    names(out) <- c('CTF_up','CTF_down','CTF_updown','TF_used')
    return(out)
}


CTF_filter <- function(CTF_result, up.fit.cutoff, down.fit.cutoff, updown.fit.cutoff, updown.tstat.cutoff){
    up <- CTF_result[[1]]
    down <- CTF_result[[2]]
    updown <- CTF_result[[3]]

    'filter by cutoff'
    up.idx <- as.numeric(up[,1])>=up.fit.cutoff
    down.idx <- as.numeric(down[,1])>=down.fit.cutoff
    updown.idx <- as.numeric(updown[,1])>=updown.fit.cutoff & as.numeric(updown[,3])>=updown.tstat.cutoff
    # print(up.idx)

    'waring message'
    if(sum(up.idx)<1){print('up.fit.cutoff is too large!')}
    if(sum(down.idx)<1){print('down.fit.cutoff is too large!')}
    if(sum(updown.idx)<1){print('updown.fit.cutoff or updown.tstat.cutoff is too large!')}

    up.filter <- as.data.frame(up[up.idx,])
    # print(up.filter)
    down.filter <- as.data.frame(down[down.idx,])
    updown.filter <- as.data.frame(updown[updown.idx,])

    if (sum(up.idx)==1){
        up.filter <- data.frame(lapply(up.filter, function(x) t(data.frame(x))))
        rownames(up.filter) <- rownames(up)[up.idx]
        colnames(up.filter) <- c('goodness_of_fit', 'Asym', 'xmid', 'scale', 'stablesite_low', 'stablesite_high')
    }
    if (sum(down.idx)==1){
        down.filter <- data.frame(lapply(down.filter, function(x) t(data.frame(x))))
        rownames(down.filter) <- rownames(down)[down.idx]
        colnames(down.filter) <- c('goodness_of_fit', 'Asym', 'xmid', 'scale', 'stablesite_low', 'stablesite_high')
    }
    if (sum(updown.idx)==1){
        updown.filter <- data.frame(lapply(updown.filter, function(x) t(data.frame(x))))
        rownames(updown.filter) <- rownames(updown)[updown.idx]
        colnames(updown.filter) <- c('goodness_of_fit_best', 'scale', 't_statistic', 'best_time_points')
    }
    # print(up.filter)
    CTF_result[[1]] = up.filter[order(as.numeric(up.filter[,1]),decreasing=T),]
    CTF_result[[2]] = down.filter[order(as.numeric(down.filter[,1]),decreasing=T),]
    CTF_result[[3]] = updown.filter[order(as.numeric(updown.filter[,1]),decreasing=T),]
    CTF_result[[1]] = CTF_result[[1]][setdiff(rownames(CTF_result[[1]]), rownames(CTF_result[[3]])), ]
    CTF_result[[2]] = CTF_result[[2]][setdiff(rownames(CTF_result[[2]]), rownames(CTF_result[[3]])), ]
    # print(CTF_result[[1]])
    return(CTF_result)
}


CTF_filter_auto <- function(CTF_result, top.pct=0.1){
    up <- CTF_result[[1]]
    down <- CTF_result[[2]]
    updown <- CTF_result[[3]]

    up.fit.cutoff <- quantile(as.numeric(up[,1]), 1-top.pct)
    down.fit.cutoff <- quantile(as.numeric(down[,1]), 1-top.pct)
    updown.fit.cutoff <- quantile(as.numeric(updown[,1]), 1-top.pct)
    updown.tstat.cutoff <- 0

    print(paste('up.fit.cutoff', up.fit.cutoff, sep=' : '))
    print(paste('down.fit.cutoff', down.fit.cutoff, sep=' : '))
    print(paste('updown.fit.cutoff', updown.fit.cutoff, sep=' : '))
    print(paste('updown.tstat.cutoff', updown.tstat.cutoff, sep=' : '))

    ret <- CTF_filter(CTF_result, up.fit.cutoff, down.fit.cutoff, updown.fit.cutoff, updown.tstat.cutoff)
    return(ret)
}

CTF_ranking <- function(expRdata, TFpath, path, up.fit.cutoff=0.5, down.fit.cutoff=0.5, updown.fit.cutoff=0.5, updown.tstat.cutoff=10){
    CTF_result <- CTF_finding(expRdata, TFpath, path) # plot fitness and tstatistic distribution
    'Core TF filter'
    CTF_result_filter <- CTF_filter(CTF_result, up.fit.cutoff, down.fit.cutoff, updown.fit.cutoff, updown.tstat.cutoff)
    return(CTF_result_filter)
}


'2. Core TF Figure'
'Plot gene expression of each cell and fitted sigmoid line'
#绘制某基因在各样本中的表达值，并用logistic function进行拟合绘制拟合曲线
gene_exp_logistic_fitness <- function(exp_list, gene, state_start, state_end, updown_best, TFtype){
    for (i in 1:length(exp_list)){rownames(exp_list[[i]]) <- toupper(rownames(exp_list[[i]]))}
    # 计算每个时间点的样本数
    num_vec <- c(0)
    num_add <- 0
    for (i in 1:length(exp_list)){
        num <- ncol(exp_list[[i]])
        num_add <- num_add+num
        num_vec <- c(num_vec,num_add)
    }
    # 某基因表达值
    exp <- c()
    color <- brewer.pal(state_end-state_start+1,'Set3')
    colors <- c()
    j=1
    for (i in state_start:state_end){
        exp <- c(exp,as.numeric(exp_list[[i]][gene,]))
        colors <- c(colors,rep(color[j],ncol(exp_list[[i]])))
        j <- j+1
    }
    exp[exp==0] <- 0.001
    data <- data.frame('x'=seq(length(exp)),'y'=exp)
    data_fit <- data[seq(num_vec[updown_best+1]),]#这个注意一下
    # 拟合
    # NM=try(nls_model <- nls(y~ SSlogis(x, Asym, xmid, scal), data_fit),silent = T) #Asym渐近线，xmid 1/2渐近线对应的x值，scale对x轴缩放
    # if(!"try-error" %in% class(NM)){
    #     goodness_of_fit <- cor(as.numeric(data_fit[,'y']),predict(nls_model)) #拟合优度
    #     coefficient <- coef(nls_model) #sigmoid曲线估计出来的3个参数
    #     pre <- predict(nls_model) #sigmoid曲线预测值
    #     plot(as.numeric(data[,'x']), as.numeric(data[,'y']), col=colors, pch=20, type = "p",cex=1.3,main=gene,ylim=c(0,max(as.numeric(data[,'y']))),xlab='sample.NO',ylab='log2(expression)')
    #     lines(as.numeric(data_fit[,'x']),pre,type='l',col=rgb(0.2,0.8,0.2,0.5), lwd=4)
    #     # smooth = smooth.spline(as.numeric(data[,'x']), as.numeric(data[,'y']), spar=1)
    #     # lines(smooth, type="l", col=rgb(0.2,0.8,0.2,0.5),lwd=4)
    # }else{
    #     print("Can not fit logistic function!")
    # }
    plot(as.numeric(data[,'x']), as.numeric(data[,'y']), col=colors, pch=20, type = "p",cex=1.3,main=gene,ylim=c(0,max(as.numeric(data[,'y']))),xlab='sample.NO',ylab='log2(expression)')
    smooth = smooth.spline(as.numeric(data[,'x']), as.numeric(data[,'y']), spar=1)
    lines(smooth, type="l", col=rgb(0.2,0.8,0.2,0.5),lwd=4)
}


'Core TF expression boxplot'
CTF.boxplot <- function(TFdata, CTF_result, savepath){
    pdf(paste(savepath, 'up.TFs.boxplot.pdf', sep='/'))
    par(mfrow=c(3,3))
    for (genei in rownames(CTF_result[[1]])){
        data <- list()
        for (i in 1:length(TFdata)){
            data[[i]] <- as.numeric(TFdata[[i]][genei,])
        }
        boxplot(data,main=genei,names=names(TFdata))
        stripchart(data,vertical=T,method='jitter',cex=0.8,pch=19,col=brewer.pal(length(names(TFdata)),'Set3'),add=T)
    }
    dev.off()

    pdf(paste(savepath, 'down.TFs.boxplot.pdf', sep='/'))
    par(mfrow=c(3,3))
    for (genei in rownames(CTF_result[[2]])){
        data <- list()
        for (i in 1:length(TFdata)){
            data[[i]] <- as.numeric(TFdata[[i]][genei,])
        }
        boxplot(data,main=genei,names=names(TFdata))
        stripchart(data,vertical=T,method='jitter',cex=0.8,pch=19,col=brewer.pal(length(names(TFdata)),'Set3'),add=T)
    }
    dev.off()

    pdf(paste(savepath, 'updown.TFs.boxplot.pdf', sep='/'))
    par(mfrow=c(3,3))
    for (genei in rownames(CTF_result[[3]])){
        # print(genei)
        data <- list()
        for (i in 1:length(TFdata)){
            data[[i]] <- as.numeric(TFdata[[i]][genei,])
        }
        boxplot(data,main=genei,names=names(TFdata),outpch=NA, las=2)
        stripchart(data,vertical=T,method='jitter',cex=0.8,pch=19,col=brewer.pal(length(names(TFdata)),'Set3'),add=T)
    }
    dev.off()   
}


'Core TF expression fitted line'
CTF.fitline.plot <- function(TFdata, CTF_result, savepath){
    if (nrow(CTF_result[[1]])>0){
        pdf(paste(savepath, 'up.TFs.fitness.pdf', sep='/'))
        par(mfrow=c(3,3))
        for (genei in rownames(CTF_result[[1]])){
            gene_exp_logistic_fitness(TFdata,genei, 1, length(TFdata), length(TFdata), 'up')
        }
        dev.off()
    }

    if (nrow(CTF_result[[2]])>0){
        pdf(paste(savepath, 'down.TFs.fitness.pdf', sep='/'))
        par(mfrow=c(3,3))
        for (genei in rownames(CTF_result[[2]])){
            gene_exp_logistic_fitness(TFdata,genei, 1, length(TFdata), length(TFdata), 'down')
        }
        dev.off()
    }

    if (nrow(CTF_result[[3]])>0){
        pdf(paste(savepath, 'updown.TFs.fitness.pdf', sep='/'))
        par(mfrow=c(3,3))
        for (genei in rownames(CTF_result[[3]])){
            times <- names(TFdata)
            for (i in 2:(length(times)-1)){
                if (CTF_result[[3]][genei,4] == times[i]){
                    gene_exp_logistic_fitness(TFdata,genei, 1, length(times), i, 'updown')
                }
            }
        }
        dev.off()
    }
}




####
####  3. Core TF-pair detection
####
'1. Find core cell fate TF-pair'
'Get TF-pair with two state'
get.up.down.pair <- function(data, genepair.result){
    upTFs <- rownames(genepair.result[[1]])
    downTFs <- rownames(genepair.result[[2]])
    alldata <- c()
    for (i in names(data)){
        alldata <- cbind(alldata, data[[i]])
    }
    rsquare <- c()
    upcol <- c()
    downcol <- c()
    coef <- c()
    for (up in upTFs){
        for (down in downTFs){
            genepair.exp <- t(alldata[c(up, down), ])
            colnames(genepair.exp) <- c('up', 'down')
            fit <- lm(down~up, data=as.data.frame(genepair.exp))
            coef <- c(coef, fit$coefficients[2])
            rsquare <- c(rsquare, summary(fit)$r.squared)
            upcol <- c(upcol, up)
            downcol <- c(downcol, down)
        }
    }
    up.down.pair <- data.frame(upgenes = upcol, downgenes = downcol, coef = coef, r2 = rsquare) # pair adjusted R square
    up.down.pair.filter <- up.down.pair[up.down.pair[,'coef']<0, ]
    up.down.pair.filter <- up.down.pair.filter[order(up.down.pair.filter[, 'r2'], decreasing=T), ]

    'plot r2 distribution(up-down)'
    plot(sort(as.numeric(up.down.pair.filter[, 'r2'])), seq(nrow(up.down.pair.filter))/nrow(up.down.pair.filter), col='#636363', pch=20, type = "l", lwd=4, 
        xlab='R2', ylab=paste('Cumulative Ratio (total gene: ', nrow(up.down.pair.filter), ')', sep=''), main='Distribution of R2(up-down)')
    return(up.down.pair.filter)
}


'Get TF-pair with one state'
get.up.updown.pair <- function(data, genepair.result){
    upTFs <- rownames(genepair.result[[1]])
    updownDF <- genepair.result[[3]]
    alldata <- c()
    time.cnt <- c()
    for (i in names(data)){
        alldata <- cbind(alldata, data[[i]])
        time.cnt <- c(time.cnt, ncol(data[[i]]))
    }
    time.cnt <- cumsum(time.cnt)

    upcol <- c()
    updowncol <- c()
    upcoef <- c()
    downcoef <- c()
    uprsquare <- c()
    downrsquare <- c()
    for (up in upTFs){
        for (i in 1:nrow(updownDF)){
            updown <- rownames(updownDF)[i]
            genepair.exp <- t(alldata[c(up, updown), ])
            colnames(genepair.exp) <- c('up', 'updown')
            updown.time <- updownDF[i, 4] # uodown 时间
            idx <- which(names(data)==updown.time)
            genepair.exp.up <- genepair.exp[seq(time.cnt[idx]), ]
            genepair.exp.down <- genepair.exp[seq(time.cnt[idx-1]+1, time.cnt[length(time.cnt)]), ]
            up.fit <- lm(updown~up, data=as.data.frame(genepair.exp.up))
            down.fit <- lm(updown~up, data=as.data.frame(genepair.exp.down))
            upcol <- c(upcol, up)
            updowncol <- c(updowncol, updown)
            upcoef <- c(upcoef, up.fit$coefficients[2])
            downcoef <- c(downcoef, down.fit$coefficients[2])
            uprsquare <- c(uprsquare, summary(up.fit)$r.squared)
            downrsquare <- c(downrsquare, summary(down.fit)$r.squared)
        }
    }
    up.updown.pair <- data.frame(upgenes = upcol, updowngenes = updowncol, upcoef = upcoef, upr2 = uprsquare , downcoef=downcoef, downr2 = downrsquare) # pair adjusted R square
    up.updown.pair <- up.updown.pair[!apply(up.updown.pair, 1, function(x) any(is.na(x))), ]
    up.updown.pair.filter <- up.updown.pair[up.updown.pair[,'upcoef']>0, ]
    up.updown.pair.filter <- up.updown.pair.filter[order(up.updown.pair.filter[,'upr2'], decreasing=T), ]

    'plot r2 distribution(up-updown)'
    plot(sort(as.numeric(up.updown.pair.filter[, 'upr2'])), seq(nrow(up.updown.pair.filter))/nrow(up.updown.pair.filter), col='#636363', pch=20, type = "l", lwd=4, 
        xlab='R2', ylab=paste('Cumulative Ratio (total gene: ', nrow(up.updown.pair.filter), ')', sep=''), main='Distribution of R2(up-updown)')
    return(up.updown.pair.filter)
}


'Merge two types cross-inhibit TF-pair'
CTFP_finding <- function(data, genepair.result, path){
    pdf(paste(path, 'R2.distribution.pdf', sep='/'))
    par(mfrow=c(2,2))
    up.down.pair <- get.up.down.pair(data, genepair.result)
    up.updown.pair <- get.up.updown.pair(data, genepair.result)
    dev.off()
    out <- list()
    up.updown.pair <- up.updown.pair[up.updown.pair[,'downcoef']<0.1, ]
    out[['CTFP_up_down']] <- up.down.pair
    out[['CTFP_up_updown']] <- up.updown.pair
    return(out)
}

CTFP_finding_type1 <- function(data, genepair.result, path){
    pdf(paste(path, 'R2.distribution.pdf', sep='/'))
    par(mfrow=c(2,2))
    up.down.pair <- get.up.down.pair(data, genepair.result)
    dev.off()
    out <- list()
    out[['CTFP_up_down']] <- up.down.pair
    return(out)
}


'Filter two types cross-inhibit TF-pair by R2 cutoff'
CTFP_filtering <- function(data, CTF, path, up.down.r2.cutoff=0.5, up.updown.r2.up.cutoff=0.1, up.updown.r2.down.cutoff='None'){
    CTFP_result <- CTFP_finding(data, CTF, path)
    up.down <- CTFP_result[[1]]
    up.updown <- CTFP_result[[2]]
    up.down.pair.filter <- up.down[up.down[,'r2']>up.down.r2.cutoff, ]
    if (up.updown.r2.down.cutoff=='None'){
        up.updown.pair.filter <- up.updown[up.updown[,'upr2']>up.updown.r2.up.cutoff, ]
        } else {
            up.updown.pair.filter <- up.updown[up.updown[,'upr2']>up.updown.r2.up.cutoff & up.updown[,'downcoef']<up.updown.r2.down.cutoff, ]
        }
    
    CTFP_result[[1]] <- up.down.pair.filter
    CTFP_result[[2]] <- up.updown.pair.filter
    return(CTFP_result)
}

CTFP_filter_type1 <- function(CTFP_result, up.down.r2.cutoff){
    up.down <- CTFP_result[[1]]
    up.down.pair.filter <- up.down[up.down[,'r2']>up.down.r2.cutoff, ]
    CTFP_result[[1]] <- up.down.pair.filter
    return(CTFP_result)
}

CTFP_filter_auto <- function(CTFP_result, top.pct=0.1){
    up.down <- CTFP_result[[1]]
    up.updown <- CTFP_result[[2]]

    up.down.r2.cutoff <- quantile(as.numeric(up.down[,'r2']), 1-top.pct)
    up.updown.r2.cutoff <- quantile(as.numeric(up.updown[,'upr2']), 1-top.pct)

    print(paste('up.down.r2.cutoff', up.down.r2.cutoff, sep=' : '))
    print(paste('up.updown.r2.cutoff', up.updown.r2.cutoff, sep=' : '))

    ret <- CTFP_filter(CTFP_result, up.down.r2.cutoff, up.updown.r2.cutoff)
    return(ret)
}


'2. Show TF-pair expression with smmoth lines'
'One TF-pair'
get.exp.data <- function(data, genepair, colors){
    alldata <- c()
    for (i in names(data)){
        alldata <- cbind(alldata, data[[i]])
    }
    genepair.exp <- t(alldata[genepair, ])
    par(mar=c(5, 5, 3, 9))
    plot(seq(nrow(genepair.exp)), genepair.exp[,1], col=colors[1], pch=20, ylab='Normalized Exp', xlab='cells', main=paste(genepair[1], genepair[2], sep='-'))
    points(seq(nrow(genepair.exp)), genepair.exp[,2], col=colors[2], pch=20)
    smoothingSpline_up = smooth.spline(seq(nrow(genepair.exp)), genepair.exp[,1], spar=1)
    smoothingSpline_updown = smooth.spline(seq(nrow(genepair.exp)), genepair.exp[,2], spar=1)
    lines(smoothingSpline_up, type="l", col=colors[1], lwd=5)
    lines(smoothingSpline_updown, type="l", col=colors[2], lwd=5)
    usr <- par("usr")
    x <- usr[2]*1.02
    y <- usr[4]*1 
    legend(x, y, legend=genepair, col=colors, pch=20, bty="o", box.lty=0, xpd=TRUE)
}


'All TF-pairs'
CTFP.fitline.plot <- function(TFdata, CTFP_result, savepath){
    up.down.pair <- CTFP_result[['CTFP_up_down']]
        if (nrow(up.down.pair)>0){
        pdf(paste(savepath, 'up.down.TFpairs.pdf', sep='/'))
        par(mfrow=c(2,1))
        for (i in 1:nrow(up.down.pair)){
            get.exp.data(data=TFdata, genepair=c(as.character(up.down.pair[i,'upgenes']), as.character(up.down.pair[i,'downgenes'])), colors=c('#ef3b2c', '#66c2a5'))
        }
        dev.off()
    }   

    up.updown.pair <- CTFP_result[['CTFP_up_updown']]
    if (nrow(up.updown.pair)>0){
        pdf(paste(savepath, 'up.updown.TFpairs.pdf', sep='/'))
        par(mfrow=c(2,1))
        for (i in 1:nrow(up.updown.pair)){
            get.exp.data(data=TFdata, genepair=c(as.character(up.updown.pair[i,'upgenes']), as.character(up.updown.pair[i,'updowngenes'])), colors=c('#ef3b2c', '#66c2a5'))
        }
        dev.off()
    }
}


'3. Order cell by ENI algorithm'
Rescale <- function (Data){
    InEC = Data
    ECNo0 = InEC[which(rowMeans(InEC) > 0), ]
    ECSC = t(apply(ECNo0, 1, scale))
    rownames(ECSC) = rownames(ECNo0)
    colnames(ECSC) = colnames(ECNo0)
    ECSC
}

PipeRCDF <- function (Data, Ndg = 20){
    t1 <- PipeR(Data, Ndg) # 每个基因基于样本的残差均值
    EC <- quantile(ecdf(t1), 1:100/100) # 累积分布
    out <- sum(EC)
}

PipeR <- function (Data, Ndg = 3, Method = "Poly"){
    V <- 1:ncol(Data) # 样本数
    if (Method == "Poly") 
        ResV <- t(sapply(1:nrow(Data), function(i) mean(residuals(lm(Data[i, 
            ] ~ poly(V, Ndg)))^2))) # 每个基因计算一个残差均值，gene*cell matrix ~ poly(cellnumber, 3)
    colnames(ResV) <- rownames(Data)
    ResV
}

ImpTC <- function (Data, Seq, condn, Ndg = 3){ #     res.imp <- ImpTC(data.sc, start.order, cond.num, Ndg = Ndg)
    if (is.null(Seq)) 
        Seq = 1:ncol(Data)
    expect_is(Seq, "integer")
    conduse <- condn[Seq]
    Od <- Seq[1:3] # 乱序后前三个样本
    Od <- Od[order(conduse[1:3])] # 对前三个样本基于时间进行排列
    condnow <- condn[Od] # 对前三个样本的时间点
    Ncol <- ncol(Data) # 样本总数
    for (i in 4:ncol(Data)) {
        tp <- Seq[i]
        tpc <- conduse[i]
        w1 <- which(condnow < tpc)
        w2 <- which(condnow > tpc)
        condleft <- ifelse(length(w1) == 0, NA, max(w1))
        condright <- ifelse(length(w2) == 0, NA, min(w2))
        if (!is.na(condleft) & !is.na(condright)) 
            mat = t(sapply(condleft:(condright - 1), function(j) c(Od[1:j], 
                tp, Od[(j + 1):(i - 1)]))) # 新的样本点及老样本点所有可能的排列，每行一种排列组合
        if (is.na(condleft)) 
            mat <- matrix(c(tp, Od[1:(i - 1)]), nrow = 1)
        if (is.na(condright)) 
            mat <- matrix(c(Od[1:(i - 1)], tp), nrow = 1)
        if (is.na(condright) & is.na(condleft)) 
            mat <- rbind(c(tp, Od[1:(i - 1)]), c(Od[1:(i - 1)], 
                tp))
        Stat = sapply(1:nrow(mat), function(j) PipeRCDF(Data[, 
            mat[j, ]], Ndg)) # 对每种顺序进行计算，计算累积错误率
        Min = which.min(Stat)
        Od = mat[Min, ]
        condnow = condn[Od]
        message("insert ", i, "th cell")
    }
    return(Od)
}


Opt2TC <- function (Data, N, Seq, condn, Ndg = 3, NCThre = 1000){ # res.out <- Opt2TC(data.sc, N, res.imp, cond.num, Ndg = Ndg, NCThre = NCThre)
    message("2-opt: ")
    Ncol <- ncol(Data)
    conduse <- condn[Seq] # 正确的顺序(离散时间)
    SeqIter <- Seq # 样本的编号
    StatIter <- PipeRCDF(Data[, Seq], Ndg) # 细胞在该排布下的累积残差
    nc = 0
    condl <- unique(conduse) # 离散时间点
    st <- sapply(1:length(condl), function(i) which(conduse == 
        condl[i])[1]) # 每个时间点的start sample
    ed <- Ncol
    if (length(condl) > 1) 
        ed <- c(sapply(1:length(condl), function(i) which(conduse == 
            condl[i])[1] - 1)[-1], Ncol) # 每个时间点end sample
    for (i in 1:N) { # 循环次数2w
        choose.l <- sample(condl, 1) # 某个时间点
        Choose <- sample(st[choose.l]:ed[choose.l], 2) # 在某个时间中选取两个样本
        Seq0 <- SeqIter # 样本编号
        Seq0[Choose[1]:Choose[2]] <- SeqIter[Choose[2]:Choose[1]]
        Stat <- PipeRCDF(Data[, Seq0], Ndg) # 随机调换两个时间点看累积残差是不是比较小，如果比较小就进行更新，防止之前陷入local minimum
        nc <- nc + 1
        if (Stat < StatIter) {
            SeqIter = Seq0
            StatIter = Stat
            nc = 0
            message("update ", i)
        }
        if (nc > NCThre) {
            message("final update:", i) # 更新1000次停止
            break
        }
    }
    Out = list(SeqIter, StatIter)
    return(Out)
}

Order_cell <- function (TFdata, CTF, savepath, Ndg = 3, N = 20000, NCThre = 1000, Seed = 1000){  # WaveCrestENI(Markers, WaveCrestExData, Conditions, N=1000) 
    Data <- c()
    for (i in names(TFdata)){
        Data <- cbind(Data, TFdata[[i]])
    }

    GeneList <- as.character(unlist(lapply(CTF, function(x) rownames(x))))

    tm.len <- as.numeric(unlist(lapply(TFdata, function(x) ncol(x))))
    CondVector <- rep(paste("t",1:length(tm.len),sep=""), tm.len)
    Conditions <- factor(CondVector, levels=paste("t",1:length(tm.len),sep=""))

    expect_is(Data, "matrix")
    set.seed(Seed)
    start.order <- sample(1:ncol(Data), ncol(Data))
    if (!is.factor(Conditions)) 
        Conditions <- factor(Conditions, levels = unique(Conditions))
    cond.num <- as.numeric(Conditions)
    if (length(setdiff(GeneList, rownames(Data))) > 0) 
        stop("some genes in GeneList are not in Data!")
    data.use <- Data[GeneList, ]
    data.sc <- Rescale(data.use) # 基于基因进行scale, 行进行scale
    res.imp <- ImpTC(data.sc, start.order, cond.num, Ndg = Ndg) # 返回的是排好序的时间点
    res.out <- Opt2TC(data.sc, N, res.imp, cond.num, Ndg = Ndg, NCThre = NCThre)
    order.Data <- Data[, res.out[[1]]]
    saveRDS(res.out[[1]], paste(savepath, 'ENIRes.rds', sep='/'))
    saveRDS(order.Data, paste(savepath, 'rank.data.rds', sep='/'))
    return(order.Data)
}


'4. Compete model for cell fate TF-pair rank'
'Calculate one gene pair fitness and create fitness figure'
compete.gene.pair <- function(alldata, gene1, gene2) {
    set.seed(123)
    # compete model for two gene
    compete <- function(time, y, parms, ...){
        with(as.list(c(parms, y)), {
            dx <- r1*x*(1-x/n1-s1*y/n2)
            dy <- r2*y*(1-y/n2-s2*x/n1)
            list(c(dx, dy))
        })
    }

    # cost function(OLS)
    cost <- function(p, data, ...) {
        yy <- p[c('x', 'y')]
        pp <- p[c('r1', 'r2', 'n1', 'n2', 's1', 's2')]
        times <- seq(nrow(data))
        out <- ode(yy, times, compete, pp)
        modCost(out, data, weight='none', ...) # try weight = "std" or "mean"
    }

    # fit compete model
    dat <- data.frame(time=seq(ncol(alldata)), x=alldata[gene1,], y=alldata[gene2, ] ) # 1.783 ######################
    times <- seq(ncol(alldata))
    errors <- c();fitnesses <- c();pvalues <- c();outs <- list();paramslist <- c()
    for (r in c(0.01, 0.1, 10)){
        r1 <- r
        r2 <- r
        n1 <- max(dat[,'x'])
        for (n2 in c(0.1, n1, n1*10, n1*100)){
            s1 <- runif(1)
            s2 <- runif(1, min=1, max=10)
            parms.init <- c(x=dat[1,'x']+0.1, y=dat[1,'y']+0.1, r1 = r1, r2 = r2, n1 = n1, n2 = n2, s1 = s1, s2 = s2)
            paramslist <- c(paramslist, parms.init)
            try_modFit <- try(fit <- modFit(f = cost, p = parms.init, data = dat, lower=rep(0,8), method='Marq'), silent = TRUE) ######################
            # print(parms.init)
            # print(class(try_modFit))
            if(!'try-error' %in% class(try_modFit)){
                    y_init <- coef(fit)[c('x', 'y')]
                    params <- coef(fit)[c('r1', 'r2', 'n1', 'n2', 's1', 's2')]
                    out <- ode(y_init, times, compete, params)
                    error <- sum((dat-out)**2)
                    fitness <- cor(c(as.numeric(out[,'x']), as.numeric(out[,'y'])), c(as.numeric(dat[,'x']), as.numeric(dat[,'y'])))
                    pvalue <- cor.test(c(as.numeric(out[,'x']), as.numeric(out[,'y'])), c(as.numeric(dat[,'x']), as.numeric(dat[,'y'])))$p.value
                    errors <- c(errors, error)
                    fitnesses <- c(fitnesses, fitness)
                    pvalues <- c(pvalues, pvalue)
                    outs <- c(outs, list(out))
            }

        }
    }
    
    # return fit results
    # print(errors)
    # print(fitnesses)

    if (!is.null(errors)){
        idx <- which(errors==min(errors, na.rm = T))[1]
        best_fitness <- fitnesses[idx]
        best_pvalue <- pvalues[idx]
        best_out <- outs[[idx]]
        # fitness
        if (best_fitness=='None'){
            best_fitness <- 'NA'
        }
        return(list(fitness = c(gene1, gene2, best_fitness, best_pvalue), data=dat, out=best_out))
    }else{
        return('NA')
    }

    # return(list(paramslist=paramslist, errors=errors, fitnesses=fitnesses))
}


# compete.gene.pair.parallel <- function(alldata, CTFP, ncore=5){
#     alldata <- alldata
#     start <- Sys.time()
#     options(warn=-1)

#     cluster <- makeCluster(ncore)
#     # print(1)
#     clusterExport(cluster, c('compete.gene.pair', 'modFit', 'ode', 'modCost'))
#     # print(2)
#     up_down <- clusterMap(cluster, compete.gene.pair, alldata, as.character(CTFP[[1]][seq(1),1]), as.character(CTFP[[1]][seq(1),2]))
#     up_updown <- clusterMap(cluster, wrapper, as.character(CTFP[[2]][,1]), as.character(CTFP[[2]][,2]))
#     # print(3)
#     # do.call('rbind', result)
#     stopCluster(cluster)

#     end <- Sys.time()
#     print(end - start)
#     options(warn=0)
#     # return(up_down)
#     return(list(up_down, up_updown))
# }


'compete TF-pair fit plot'
compete.fit.plot <- function(fit_ret, savepath, pdfnm){
    pdf(paste(savepath, pdfnm, sep='/'))
    par(mfrow=c(2,2))
    for (i in 1:length(fit_ret)){
        # plot
        if (length(fit_ret[[i]][[1]])>1){
            gene1 <- fit_ret[[i]][[1]][1]
            gene2 <- fit_ret[[i]][[1]][2]
            dat <- fit_ret[[i]][[2]]
            best_out <- fit_ret[[i]][[3]]
            par(mar=c(5, 5, 3, 9))
            plot(dat[,'time'], dat[,'x'], col='#fc9272', pch=16, xlab='Time', ylab='Expression', main=paste(paste(gene1, gene2, sep='-'), ' Compete Model'))
            lines(best_out[,'time'], best_out[,'x'], col='#de2d26', lwd=4)
            points(dat[,'time'], dat[,'y'], col='#a1d99b', pch=16)
            lines(best_out[,'time'], best_out[,'y'], col='#31a354', lwd=4)
            usr <- par("usr")
            x <- usr[2]*1.02
            y <- usr[4]*1 
            legend(x, y, c(gene1, gene2), pch=c(16,16), col=c('#fc9272', '#a1d99b'), box.lty=0, cex=1, xpd=TRUE)
        }
      
    }
    dev.off()
}


'Get TF-pair fitness and plot fitness line'
CTFP_ranking <- function(CTFP, compete_ret, savepath){
    # fit compete model
    # compete_ret <- compete.gene.pair.parallel(alldata=alldata, CTFP=CTFP, ncore=ncore)
    up_down <- compete_ret[[1]]
    up_updown <- compete_ret[[2]]
    # plot fit line
    compete.fit.plot(fit_ret=up_down, savepath=savepath, pdfnm='compete.TFpair.fitness.up_down.pdf')
    compete.fit.plot(fit_ret=up_updown, savepath=savepath, pdfnm='compete.TFpair.fitness.up_updown.pdf')

    get.fitness.pvalue <- function(ret){
        fitness <- c();pvalue <- c()
        for (i in 1:length(ret)){
            if (length(ret[[i]][[1]])>1){
                fitness <- c(fitness, as.numeric(ret[[i]][[1]][3]))
                pvalue <- c(pvalue, as.numeric(ret[[i]][[1]][4]))
            }else{
                fitness <- c(fitness, 0)
                pvalue <- c(pvalue, 1)
            }
        }
        return(list(fitness, pvalue))
    }
    updown_sta <- get.fitness.pvalue(up_down)
    upupdown_sta <- get.fitness.pvalue(up_updown)
    updown_fit <- updown_sta[[1]]
    upupdown_fit <- upupdown_sta[[1]]
    updown_pval <- updown_sta[[2]]
    upupdown_pval <- upupdown_sta[[2]]
    # add TF-pair fitness
    CTFP$CTFP_up_down <- CTFP$CTFP_up_down[,1:2]
    CTFP$CTFP_up_updown <- CTFP$CTFP_up_updown[,1:2]
    colnames(CTFP$CTFP_up_down) <- c('TCFM', 'SCFM')
    colnames(CTFP$CTFP_up_updown) <- c('TCFM', 'ACFM')
    CTFP$CTFP_up_down <- as.data.frame(cbind(CTFP$CTFP_up_down, Fitness = updown_fit, Pval = updown_pval))
    CTFP$CTFP_up_updown <- as.data.frame(cbind(CTFP$CTFP_up_updown, Fitness = upupdown_fit, Pval = upupdown_pval))
    CTFP$CTFP_up_down <- CTFP$CTFP_up_down[order(CTFP$CTFP_up_down[,'Fitness'], decreasing=TRUE),]
    CTFP$CTFP_up_updown <- CTFP$CTFP_up_updown[order(CTFP$CTFP_up_updown[,'Fitness'], decreasing=TRUE),]
    rownames(CTFP$CTFP_up_down) <- seq(nrow(CTFP$CTFP_up_down))
    rownames(CTFP$CTFP_up_updown) <- seq(nrow(CTFP$CTFP_up_updown))
    # return CTFP
    return(CTFP)
}









# 'Calculate all gene pairs fitness and fitness figure'
# compete.gene.pairs <- function(TFdata, gene_pairs, savepath, pdfnm) {
#     alldata <- c()
#     for (i in names(TFdata)){
#         alldata <- cbind(alldata, TFdata[[i]])
#     }

#     options(warn=-1)
#     gene_pair_fitness <- c()
#     pdf(paste(savepath, paste('Compete.Model.', pdfnm, '.pdf', sep=''), sep='/'))
#     par(mfrow=c(2,2))
#     for (i in 1:nrow(gene_pairs)){
#     # for (i in 1:2){
#         ret <- compete.gene.pair(alldata, as.character(gene_pairs[,1][i]), as.character(gene_pairs[,2][i]))
#         gene_pair_fitness <- rbind(gene_pair_fitness, ret)
#     }
#     dev.off()

#     gene_pair_fitness <- as.data.frame(gene_pair_fitness)
#     rownames(gene_pair_fitness) <- seq(nrow(gene_pair_fitness))
#     colnames(gene_pair_fitness) <- c('TF1', 'TF2', 'fitness')
#     options(warn=0)
#     return(gene_pair_fitness)
# }


# 'Rank TF-pair fitness of compete model'
# CTFP_rank <- function(TFdata, CTFP_result_filter, savepath){
#     gene_pairs_up_down <- CTFP_result_filter$CTFP_up_down
#     gene_pairs_up_updown <- CTFP_result_filter$CTFP_up_updown
#     up_down_fitness <- compete.gene.pairs(TFdata, gene_pairs_up_down, savepath, 'up_down')
#     up_updown_fitness <- compete.gene.pairs(TFdata, gene_pairs_up_updown, savepath, 'up_updown')
#     CTFP_result_filter$CTFP_up_down <- as.data.frame(cbind(as.matrix(gene_pairs_up_down), fitness = as.numeric(as.character(up_down_fitness[, 'fitness']))))
#     CTFP_result_filter$CTFP_up_updown <- as.data.frame(cbind(as.matrix(gene_pairs_up_updown), fitness = as.numeric(as.character(up_updown_fitness[, 'fitness']))))
#     CTFP_result_filter$CTFP_up_down <- CTFP_result_filter$CTFP_up_down[order(CTFP_result_filter$CTFP_up_down[,'fitness'], decreasing=TRUE),]
#     CTFP_result_filter$CTFP_up_updown <- CTFP_result_filter$CTFP_up_updown[order(CTFP_result_filter$CTFP_up_updown[,'fitness'], decreasing=TRUE),]
#     return(CTFP_result_filter)
# }


# 'Permutation of genelist1 and genelist2'
# gene.perm <- function(genevec1, genevec2){
#     gene_pairs <- c()
#     for (gene1 in genevec1){
#         for (gene2 in genevec2){
#             gene_pairs <- rbind(gene_pairs, c(gene1, gene2))
#         }
#     }
#     gene_pairs <- as.data.frame(gene_pairs)
#     return(gene_pairs)
# }

















    

