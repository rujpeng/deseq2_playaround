## load Annotations
library("GenomicFeatures")
txdb = makeTxDbFromGFF("genome.gtf",format="gtf")
ebg = exonsBy(txdb,by="gene") ## this is to group exons by gene


## load read counts from bam files
## here I compare only control and knockout
library("Rsamtools")
bam_col_1 = BamFile("Col_1.bam") ## control
bam_col_2 = BamFile("Col_2.bam") ## control
bam_col_3 = BamFile("Col_3.bam") ## control
bam_ox_1  = BamFile("OX_1.bam") ## knockout
bam_ox_2  = BamFile("OX_2.bam") ## knockout
bam_ox_3  = BamFile("OX_3.bam") ## knockout
bamfiles  = BamFileList (bam_col_1,bam_col_2,bam_col_3,bam_ox_1,bam_ox_2,bam_ox_3)

## experimental design
expdesign = read.csv("knockout_design.csv",row.names=1,sep=",")

## count reads using summarizeOverlaps
library("GenomicAlignments")
se = summarizeOverlaps(features=ebg,reads=bamfiles, ignore.strand=TRUE)
counts = assay(se)
countgenenames   = gsub("[.][1234567890]", "", row.names(counts))
rownames(counts) = countgenenames

## filtering counts, remove extremely low-read genes
medianCountByGroup = t(apply(counts,1,tapply,expdesign,median))
maxMedian = apply(medianCountByGroup,1,max)
counts_filtered = counts[maxMedian >= 0,] ## may be no filtering first?

## now DESeq2
## first normalize between samples
library("DESeq2")
cds = DESeqDataSetFromMatrix(countData=counts_filtered,colData=expdesign,
                             design= ~ condition)
cds = estimateSizeFactors(cds) ## normalization between samples
cds = estimateDispersions(cds)

pdf(file="dispersion_knockout.pdf")
plotDispEsts(cds)
dev.off()

## then report differentially expressed genes
cds = DESeq(cds)
res = results(cds)
num = sum(res$padj < 0.05, na.rm=T)
print(num)

resSigind = res[ which(res$padj < 0.05 & res$log2FoldChange > 0), ]
resSigrep = res[ which(res$padj < 0.05 & res$log2FoldChange < 0), ]
resSig = rbind(resSigind, resSigrep)

print(rownames(resSigind))
