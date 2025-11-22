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
# DESeqDataSet object for Txi Object
```{R}
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(
    txi,colData = samples,design = ~ condition
)
```
# create coldata
```{R}
coldata <- samples
coldata$files <- files
coldata$names <- coldata$run
```

# create DESeq object
```{R}
library("tximeta")
se <- tximeta(coldata)
ddsTxi_alt <- DESeqDataSet(se, design = ~ condition)
```

# remove low expressed genes, keep highly expressed genes
```{R}
dds <- ddsTxi_alt
my_rule <- 3
keep <- rowSums(counts(dds) >= 10) >= my_rule
keep 
table(keep)
dds_keep <- dds[keep,]
dds_de <- DESeq(dds_keep)
res_de <- results(dds_de)

# S4 object to data frame
res_de_df <- data.frame(res_de)
# conda install conda-forge::r-devtools
# R environment activate by typing R and Entre
# devtools::install_github("BioSenior/ggVolcano")

# https://github.com/BioSenior/ggVolcano

library(ggVolcano)

data_df <- add_regulate(
    res_de_df, log2FC_name = "log2FoldChange",
    fdr_name = "padj",log2FC = 1, fdr = 0.05
)

data_df$row <- rownames(data_df)
ggvolcano(data_df,
    x = "log2FoldChange", y = "padj",
    label = "row", label_number = 10, 
    output = FALSE)
```
