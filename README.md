# Salmon-tximport-DESEq2-Pipeline-Ver.1
Custom R and Salmon pipeline for differential gene expression analysis in my Senior honor thesis.

# PART 1: Install R, R-Studio, Linux/WSL, Ubuntu
1. Install R (Cran) - Choose Duke University, Durham NC Mirror
2. Install R Studio IDE (Posit)
3. Install WSL (for Windows 10 ver.2004 or later & Windows 11) OR use Linux:
Open PowerShell and run as administrator.
      Type:
```
wsl --install
```

5. Ubuntu (opened in WSL or can be opened as its own app) will ask: “Enter new UNIX username.”; Choose any username.
6. Then choose a password; When typing the password, nothing appears on the screen; this is normal in Linux. Type it and press Enter.

At this point, you should see a prompt like: “alex@DESKTOP:~$”

7. Type:
```
sudo apt update
```

9. Enter your newly created password.
10. This updates your Linux packages;Type:
```
sudo apt upgrade -y
```
    
# PART 2: Install Salmon inside Ubuntu  
From here you should be typing all these codes after a prompt that looks similar to “alex@DESKTOP:~$” ← this signals that you are using Ubuntu 

1. Install Salmon
 ```
sudo apt install salmon
 ```
2. Verify that salmon is installed
 ```
salmon --version
 ```
# PART 3: Using Salmon inside Ubuntu 
From here you should be typing all these codes after a prompt that looks similar to “alex@DESKTOP:~$” ← this signals that you are using Ubuntu 

1. Change the path to the folder you want to save all documents. 
Your prompt should look something like → alex@DESKTOP: /mnt/c/Users/Lab/Downloads/WSL$

2. Type:
```
cd /mnt/c/Users/Lab/Downloads/WSL
```

3. Check if your FASTQ files are located in the path you created above (and therefore the computer can access it): 
You should have something like the following in your folder: (sample1_R1.fastq.gz, sample1_R2.fastq.gz, sample2.fastq.gz)
Sample 1 is a multi-end file and sample 2 is a single-end file; I believe that we only have single-end files.

4. Type:
```
ls *.fastq.gz
```

5. Download a transcriptome FASTA and GTF file from Ensembl (https://www.ensembl.org/Homo_sapiens/Info/Index)
   FASTA and GTF files MUST be from the same version

This downloads the FASTA file into the folder you pathed above:
```
wget ftp://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
    gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
```
    
This downloads the GTF file into the folder you pathed above. This will be used later in R Studio: 
```
wget https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz
    gunzip Homo_sapiens.GRCh38.116.gtf.gz
```

7. Create an Index; This creates a folder: salmon_index
 ```
salmon index \
 -t Homo_sapiens.GRCh38.cdna.all.fa \
 -i salmon_index \
 -k 31
 ```

7. Single-end FASTQ command; Don't forget to change the title of FASTQ files in the following code to your sample name; Type: 
 ```
salmon quant \
-i salmon_index \
-l A \
-r CJQPG3_8_8_shGPR116-SS-1.fastq.gz \
-p 8 \
-o CJQPG3_8_8_shGPR116-SS-1_quant
 ``` 
8. Repeat for each sample
# Each output folder will contain a quant.sf THIS IS AN IMPORTANT FILE



# PART 4: Using RStudio: Make Gene Matrix Count
If you have already installed the packages, you do not need to reinstall.
You MUST still type the libraries into the console and reinstall GTF file if you restart RStudio.

1. Open RStudio
2. Install packages; Type:
```
install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("tximportData")
install.packages("jsonlite")
```

3. Check if you have tximport, readr, and Annotation Dbi installed; Type:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("tximport", "readr", "AnnotationDbi"))
Check if you have rtracklayer package installed; Type: 
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("rtracklayer")
```

# From here, you will need to redo the following if you exit out of RStudio.

4. Import the GTF file into Rstudio; Change the path to where your GTF file is located; Type: 
```
library(rtracklayer) 
gtf <- import("/path/to/this/document/Homo_sapiens.GRCh38.116.gtf") 
```

5. Import Salmon results into RStudio; Change path to desired sample; Type:
```
library(tximport)
files <- file.path("C:/Users/… path to sample 1’s quant.sf file", "quant.sf")
names(files) <- "sample1"
txi <- tximport(files, type="salmon", txOut=TRUE)
```

6. Type:
```
tx2gene <- unique(data.frame(
transcript_id = mcols(gtf)$transcript_id,
gene_id = mcols(gtf)$gene_id
))

tx2gene <- na.omit(tx2gene)
head(tx2gene) 
```

7. Type: 
```
library(tximport)
files <- c(sample1 = "path/to/quant.sf")
```

8. Type:
```
txi <- tximport(
files,
type = "salmon",
tx2gene = tx2gene,
ignoreTxVersion = TRUE
)
```

9. This creates the gene count matrix; Type:
```
counts <- txi$counts
```

10. This downloads the gene count matrix as a .csv file into your desired folder; Type:
```
setwd("C:/what/ever/folder/path/you/want/to/upload/this/file")
write.csv(counts, "gene_counts.csv")
```

11. Merge all csv files into one complete file; Change the file location (line 177) and desired location of csv file (line 197) Type:

```
library(readr)
library(dplyr)
library(purrr)

files <- list.files(
    "C:/Users/Lab/Downloads/WSL/Counts.csv-files",
    pattern = "\\.csv$",
    full.names = TRUE
)

count_list <- lapply(files, function(file){
    df <- read_csv(file)
    names(df)[1] <- "GeneID"
    sample_name <- tools::file_path_sans_ext(basename(file))
    names(df)[2] <- sample_name
    df
})

count_matrix <- reduce(count_list, full_join, by = "GeneID")

count_matrix[is.na(count_matrix)] <- 0

write.csv(
    count_matrix,
    "C:/Users/YourName/Documents/RNAseq/combined_counts.csv",
    row.names = FALSE)
```

# PART 5: Using RStudio: DESeq2 
If you have already installed the packages, you do not need to reinstall; you MUST still type the libraries into the console.
    
1. Load Libraries required to run differential sequencing and making volcano plots; change path to combined counts (line 211); Type:
```
install.packages("BiocManager")
BiocManager::install("EnhancedVolcano") 
library(DESeq2)
library(EnhancedVolcano)
Load your combined_counts.csv into RStudio
counts <- read.csv(
 "C:/path/to/your/combined_counts.csv",  
    row.names = 1,
    check.names = FALSE)
```

2. Create the sample metadata 
This assigns the numbers with its experimental condition; therefore allowing the computer to know what gene counts are for what experimental condition;
Edit for your experimental condition; Type

```
sample_names <- colnames(counts)

sampleInfo <- data.frame(
row.names = sample_names,

knockdown = ifelse(
grepl("shSCR", sample_names),
"shSCR",
"shGPR116"),

condition = ifelse(
grepl("-static", sample_names),
"static",
"SS")) 
```


3. Add the batch differences to the sampleInfo metadata;
This will add a new column to the table above that labels each sample with the batch it belongs to (CJQPG3 vs. X6CM82); Type:

```
sampleInfo$batch <- ifelse(
grepl("^CJQPG3", rownames(sampleInfo)),
"CJQPG3",
"X6CM82")

sampleInfo$batch <- factor(sampleInfo$batch) 
```

4. To verify the sample metadata:
Type:
```
sampleInfo
```
Type: 
```
colnames(counts)
```

# 5. Make factors (IMPORTANT for DESeq2); Type:
```
sampleInfo$knockdown <- factor(sampleInfo$knockdown, levels = c("shSCR","shGPR116"))
sampleInfo$condition <- factor(sampleInfo$condition, levels = c("static","SS"))
```

6.Build DESeq2 dataset (interaction model); Type:

```
dds <- DESeqDataSetFromMatrix(
countData = round(counts),
colData = sampleInfo,
design = ~ batch + condition + knockdown + condition:knockdown)
```

7. Remove very-low count genes; Type:
```
dds <- dds[rowSums(counts(dds)) >= 10, ]
```
8. Run DESeq2; Type:
```
dds <- DESeq(dds)
```

9. Check coefficient names (IMPORTANT); these names will be used for later codes; Type:
```
resultsNames(dds)
```

10. Create table for the results of DESeq2 of STATIC: shGPR116 vs shSCR; Type:
```
res_static <- results(dds, name = "knockdown_shGPR116_vs_shSCR")
```

11. Create a table for the results of DESeq2 of SS: shGPR116 vs shSCR; Type: 
```
res_ss <- results(dds, contrast = list(c("knockdown_shGPR116_vs_shSCR", "conditionSS.knockdownshGPR116")))
```

12. Log fold shrinkage for static and SS volcano plot; Type:

```
BiocManager::install("ashr") 
library(ashr) 

res_static <- lfcShrink( dds, coef = "knockdown_shGPR116_vs_shSCR", type = "ashr" ) 
res_ss <- lfcShrink( dds, contrast = list(c("knockdown_shGPR116_vs_shSCR", "conditionSS.knockdownshGPR116")), type = "ashr" ) 
```

13. Save the tables of results from DESeq2 as .csv files; Type:

```
write.csv(as.data.frame(res_static), "DEGs_static.csv")
write.csv(as.data.frame(res_ss), "DEGs_SS.csv") 
```

# PART 6: Using RStudio: Creating Volcano Plots

1. Convert all ENSG codes to the conventional names for the volcano plots; (ex. ENSG00000069122 → ADGRF5); Type:
Install the library first:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db") 
```

2. This now then changes the ENSG codes to the conventional names
This takes the DESeq2 object (res_static and res_ss) and convert its row_names (which are all ENSG gene IDs) into its “symbol” aka conventional name	

```
library(org.Hs.eg.db)
library(AnnotationDbi)

res_static$symbol <- mapIds(
org.Hs.eg.db,
keys = rownames(res_static),
column = "SYMBOL",
keytype = "ENSEMBL",
multiVals = "first") 

res_ss$symbol <- mapIds(
org.Hs.eg.db,
keys = rownames(res_ss),
column = "SYMBOL",
keytype = "ENSEMBL",
multiVals = "first"
) 
```

3. This step is for if you want to check ENSG codes for their conventional names; Type:

```
library(org.Hs.eg.db)
library(AnnotationDbi)
 
select(org.Hs.eg.db, keys = "GPR116", columns = c("ENSEMBL", "SYMBOL", "ALIAS"), keytype = "ALIAS")
```


4. Create Volcano plot with EnhancedVolcano function; Example Type: 
```
EnhancedVolcano(
as.data.frame(res_static),
lab = res_static$symbol,
x = "log2FoldChange",
y = "padj",
title = "shGPR116 vs shSCR (Static)",
pCutoff = 0.05,
FCcutoff = 1) 
```
5. Repeat with other conditions



# PART 7: Using RStudio: Creating .csv files for MOST upregulated/downregulated from volcano plots

1. Convert to data frame; Type:
```
“res_static_df” / “res_ss_df” are the names for these new data frames
“as.data.frame(... “ is a function that converts the DESeq2 objects (res_static and res_ss) into a data frame that can be converted into a .csv file

res_static_df <- as.data.frame(res_static)
res_ss_df <- as.data.frame(res_ss) 
```

2. Create a list of the most upregulated/downregulated genes of each volcano plot; Type:

Upregulated in static (Example):
```
up_static <- subset(
    res_static_df,
    padj < 0.05 &
    log2FoldChange > 1)
```

3. Repeat for other comparisons


4. This code makes the most significant genes (smallest adjusted p-value) appear first; Type:
```
up_static <- up_static[order(up_static$padj), ]
down_static <- down_static[order(down_static$padj), ]

up_ss <- up_ss[order(up_ss$padj), ]
down_ss <- down_ss[order(down_ss$padj), ]
```

5. Save in csv file; Type:
```
write.csv(up_static,
"Upregulated_Static.csv",
row.names = FALSE)
```
6. Repeat for all csv made
