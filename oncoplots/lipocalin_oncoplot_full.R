# libraries
library(maftools)
library(plyr)
library(tidyverse)

#todo create lipocalin gene list
lipocalin_genes <- c("OBP2A", "LCN1", "LCN6", "LCN8", "LCN9", "LCN10", "LCN12", "LCN15", "LCNL1", "PTGDS",
                     "OBP2B", "LCN1P1", "LCN2", "AMBP", "ORM1", "ORM2", "C8G", "PAEP", 
                     "ApoD", "ApoM", "CRABP1", "CRBP2", 
                     "FABP1", "FABP12", "FABP2",  "FABP3", 
                     "FABP4",  "FABP5", "FABP6", "FABP7", 
                     "FABP9", "RBP1", "RBP2", "RBP4", "RBP5", 
                     "RBP7",  "PMP2")

# read snv and csv data
oncoprint_data <- read.csv("../../../datasets/oncoprint_data_snv.csv", header = T)

# Rename the column 'hugoGeneSymbol' to 'Hugo_Symbol'
colnames(oncoprint_data)[colnames(oncoprint_data) == "hugoGeneSymbol"] <- "Hugo_Symbol"
colnames(oncoprint_data)[colnames(oncoprint_data) == "chr"] <- "Chromosome"
colnames(oncoprint_data)[colnames(oncoprint_data) == "startPosition"] <- "Start_Position"
colnames(oncoprint_data)[colnames(oncoprint_data) == "endPosition"] <- "End_Position"
colnames(oncoprint_data)[colnames(oncoprint_data) == "referenceAllele"] <- "Reference_Allele"
#colnames(oncoprint_data)[colnames(oncoprint_data) == "patientId"] <- "Tumor_Seq_Allele2"
# you can use the variant Allele as tumor_seq_allele2
colnames(oncoprint_data)[colnames(oncoprint_data) == "variantAllele"] <- "Tumor_Seq_Allele2"
colnames(oncoprint_data)[colnames(oncoprint_data) == "mutationType"] <- "Variant_Classification"
#colnames(oncoprint_data)[colnames(oncoprint_data) == "variantAllele"] <- "Variant_Type"
# you can use alterationTypee for variant_Type
colnames(oncoprint_data)[colnames(oncoprint_data) == "alterationType"] <- "Variant_Type"
colnames(oncoprint_data)[colnames(oncoprint_data) == "sampleId"] <- "Tumor_Sample_Barcode"

# replace MUTATION_EXTENDED in variant type with SNV
oncoprint_data[oncoprint_data$Variant_Type == "MUTATION_EXTENDED",]$Variant_Type <- "SNP"
#oncoprint_data[oncoprint_data$Variant_Type == "",]$Variant_Type <- "SV"

scna_data <- read.csv("../../../datasets/cnv_data/cna_obp_data.csv", header = T)
scna_data$samples <- str_replace_all(scna_data$samples, "\\.", "-")

oncoprint_data_ov <- oncoprint_data[oncoprint_data$Tumor_Sample_Barcode %in% scna_data$samples,] 

scna_data$cnv <- mapvalues(scna_data$cnv, c(-2, -1, 0, 1, 2), c("DeepDel", "Del",
                                                              "wt", "ShallowAmp",
                                                              "Amp"))
scna_data <- scna_data[scna_data$cnv != "wt",]
scna_data <- scna_data[c("samples", "Hugo_Symbol", "tumor_type", "cnv")]
colnames(scna_data) <- c("Tumor_Sample_Barcode", "Gene", "Tumor_type", "CN")
scna_data <- scna_data[c("Gene", "Tumor_Sample_Barcode", "CN", "Tumor_type")]

# filter scna and snv data to only contain lipocalins
scna_data <- scna_data[tolower(scna_data$Gene) %in% tolower(lipocalin_genes),]
oncoprint_data_ov <- oncoprint_data_ov[tolower(oncoprint_data_ov$Hugo_Symbol) %in% tolower(lipocalin_genes),]

# filter Del and shallowAmp from scna data
scna_data_filt <- scna_data[!(scna_data$CN %in% c("Del", "ShallowAmp")),]

clin_data = scna_data[c("Tumor_Sample_Barcode", "Tumor_type")]

clin_data = clin_data[!duplicated(clin_data),]

maf_object <- read.maf(maf = oncoprint_data, cnTable = scna_data, clinicalData = clin_data)
maf_object_filt <- read.maf(maf = oncoprint_data, cnTable = scna_data_filt, clinicalData = clin_data)

vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Set1')

names(vc_cols) = c("ShallowAmp", "Del", "Amp", "DeepDel", "Missense", "Truncating", "Translation_Start_Site",
                   "In-frame")

pdf("oncoprint_lipocalins_full.pdf", width = 15)
oncoplot(maf_object, includeColBarCN = TRUE, colors = vc_cols, top=50, bgCol = "white", 
         cBioPortal = TRUE, clinicalFeatures = "Tumor_type")
dev.off()

pdf("oncoprint_lipocalins_full.png")
oncoplot(maf_object, includeColBarCN = TRUE, colors = vc_cols, top=50, bgCol = "white",
         cBioPortal = TRUE)
dev.off()

names(vc_cols) = c("Amp", "DeepDel", "Missense", "Truncating", "Translation_Start_Site",
                   "In-frame", "Multi_Hit")

pdf("oncoprint_lipocalins_deep_events.pdf", width = 12)
oncoplot(maf_object_filt, includeColBarCN = TRUE, colors = vc_cols, top=50, bgCol = "white",
         cBioPortal = TRUE, clinicalFeatures = "Tumor_type")
dev.off()

