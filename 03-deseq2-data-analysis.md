# activate the conda environment
```{CMD}
conda activate r4.5
```
# start R in r4.5
```{CMD}
R
```

# inside R
```{R}
library("DESeq2")
library("DESeq2)
```


# DESeq2 Tutorial
```{R}
library("tximport")
library("readr")
library("tximportData")
# locating the directory or folder name extdata
dir <- system.file("extdata", package="tximportData")
list.files(path=dir)
# Metadata
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)

rownames(samples)
colnames(samples)

rownames(samples) = samples$run
# rownames(samples) <- samples$run
# samples$run -> rownames(samples)
# a = 2

samples$condition <- factor(rep(c("A","B"),each=3))
rownames(samples) <- samples$run
samples[,c("pop","center","run","condition")]
```
# loading salmon quantified files (quant)
```{R}
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- samples$run

# transcript to gene map
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))

#loading salmon quant files
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
```
