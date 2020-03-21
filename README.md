# CTFPdetector Overview
CTFPdetector is an R script (R version > 3.6.1) for cell fate transcription factor pairs detection in lineage commitment process. CTFPdetector aims at the cell directional differentiation system. It focuses on cells in transition states, using single-cell transcriptome data identifies mutually inhibiting transcription factor pairs which determine different cell fates.

# Installation
> 1. Download R scripts from github
> 2. Source the script

```R
source('CFTF.function.r')
```

# Quick Start
**Input file**: A count matrix of single-cell RNA-seq. The first row is time point imformation, others are genes. Columns are cells.

**Output file**: A TF pairs list. It contains transcription factor pairs, fitness of competition model and the p-value of the fitness.

1. Data preprocessing
```R
# Split count matrix by time points
> rawdata <- split.data.bytime(rawdata)
# Cell filtering
> rawdata <- filter.cells.by.MT.pect(rawdata, SAVE.path)
# Data normalization
> rawdata <- get.norm.data(rawdata, SAVE.path)
# Get more variable TFs
> TFdata <- get.variable.TF(rawdata, TFlist, SAVE.path, min.cutoff=0.1)
```

2. Core TF marker filtering
```R
> sim.TF.data <- CTF_filtering(sim.TF.data)
```

3. Core TF marker ranking
```R
> CTF_result_filter <- CTF_ranking(expRdata=sim.TF.data, TFpath=TF.path, path=SAVE.path)
```

4. Recover cell order by ENI algorithm
```R
# Cell order by ENI(extended nearest-insertion) algorithm
> rank.data <- Order_cell(TFdata=sim.TF.data, CTF=CTF_result_filter, savepath=SAVE.path, N=5000, NCThre=500)
```

5. Core TF pairs detection
```R
# Core TF pairs filtering
> CTFP_result_filter <- CTFP_filtering(data=sim.TF.data, CTF=CTF_result_filter, path=SAVE.path)

# Multithreading for speed up
> options(warn=-1)
> cluster <- makeCluster(20)
> clusterExport(cluster, c('rank.data', 'compete.gene.pair', 'modFit', 'ode', 'modCost'))
> wrapper <- function(gene1, gene2) {
>     try(compete.gene.pair(rank.data, gene1, gene2))
> }
> result1 <- clusterMap(cluster, wrapper, as.character(CTFP_result_filter[[1]][,1]), as.character(CTFP_result_filter[[1]][,2]))
> result2 <- clusterMap(cluster, wrapper, as.character(CTFP_result_filter[[2]][,1]), as.character(CTFP_result_filter[[2]][,2]))
> compete_ret <- list(result1, result2)
> stopCluster(cluster)
> options(warn=0)

# Rank TF pairs on competition model fitness
> CTFP_final <- CTFP_ranking(CTFP=CTFP_result_filter, compete_ret=compete_ret, savepath=SAVE.path)
```















