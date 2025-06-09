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

scna_data <- read.csv("../../../datasets/cnv_data/ov_cna_obp_data.csv", header = T)
scna_data$X <- str_replace_all(scna_data$X, "\\.", "-")

oncoprint_data_ov <- oncoprint_data[oncoprint_data$Tumor_Sample_Barcode %in% scna_data$X,] 

scna_data <- scna_data %>% pivot_longer(names_to = "Gene", values_to = "CN", cols = -X)

scna_data$CN <- mapvalues(scna_data$CN, c(-2, -1, 0, 1, 2), c("DeepDel", "Del",
                                                        "wt", "ShallowAmp",
                                                        "Amp"))
scna_data <- scna_data[scna_data$CN != "wt",]
colnames(scna_data) <- c("Sample_name", "Gene", "CN")
scna_data <- scna_data[c("Gene", "Sample_name", "CN")]

##scna_data <- scna_data[scna_data$Gene %in% unique(oncoprint_data$Hugo_Symbol),]

# filter scna and snv data to only contain lipocalins
scna_data <- scna_data[tolower(scna_data$Gene) %in% tolower(lipocalin_genes),]
oncoprint_data_ov <- oncoprint_data_ov[tolower(oncoprint_data_ov$Hugo_Symbol) %in% tolower(lipocalin_genes),]

# filter Del and shallowAmp from scna data
scna_data_filt <- scna_data[!(scna_data$CN %in% c("Del", "ShallowAmp")),]

maf_object <- read.maf(maf = oncoprint_data, cnTable = scna_data)
maf_object_filt <- read.maf(maf = oncoprint_data, cnTable = scna_data_filt)

vc_cols = RColorBrewer::brewer.pal(n = 12, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'DeepDel',
  'Del',
  'ShallowAmp',
  'Amp'
)

pdf("oncoprint_lipocalins_full_ov.pdf")
oncoplot(maf_object, includeColBarCN = TRUE, colors = vc_cols, top=50)
dev.off()

pdf("oncoprint_lipocalins_deep_events_ov.pdf")
oncoplot(maf_object_filt, includeColBarCN = TRUE, colors = vc_cols, top=50)
dev.off()

scna_data <- read.csv("../../../datasets/cnv_data/breast_cna_obp_data.csv", header = T)
scna_data$X <- str_replace_all(scna_data$X, "\\.", "-")

oncoprint_data_breast <- oncoprint_data[oncoprint_data$Tumor_Sample_Barcode %in% scna_data$X,] 

scna_data <- scna_data %>% pivot_longer(names_to = "Gene", values_to = "CN", cols = -X)

scna_data$CN <- mapvalues(scna_data$CN, c(-2, -1, 0, 1, 2), c("DeepDel", "Del",
                                                              "wt", "ShallowAmp",
                                                              "Amp"))
scna_data <- scna_data[scna_data$CN != "wt",]
colnames(scna_data) <- c("Sample_name", "Gene", "CN")
scna_data <- scna_data[c("Gene", "Sample_name", "CN")]

#scna_data <- scna_data[scna_data$Gene %in% unique(oncoprint_data$Hugo_Symbol),]

# filter scna and snv data to only contain lipocalins
scna_data <- scna_data[tolower(scna_data$Gene) %in% tolower(lipocalin_genes),]
oncoprint_data_breast <- oncoprint_data_breast[tolower(oncoprint_data_breast$Hugo_Symbol) %in% tolower(lipocalin_genes),]

# filter Del and shallowAmp from scna data
scna_data_filt <- scna_data[!(scna_data$CN %in% c("Del", "ShallowAmp")),]

maf_object <- read.maf(maf = oncoprint_data, cnTable = scna_data)
maf_object_filt <- read.maf(maf = oncoprint_data, cnTable = scna_data_filt)

pdf("oncoprint_lipocalins_full_breast.pdf")
oncoplot(maf_object, includeColBarCN = TRUE, colors = vc_cols, top=50)
dev.off()

pdf("oncoprint_lipocalins_deep_events_breast.pdf")
oncoplot(maf_object_filt, includeColBarCN = TRUE, colors = vc_cols, top=50)
dev.off()

scna_data <- read.csv("../../../datasets/cnv_data/melanoma_cna_obp_data.csv", header = T)
scna_data$X <- str_replace_all(scna_data$X, "\\.", "-")

oncoprint_data_melanoma <- oncoprint_data[oncoprint_data$Tumor_Sample_Barcode %in% scna_data$X,] 

scna_data <- scna_data %>% pivot_longer(names_to = "Gene", values_to = "CN", cols = -X)

scna_data$CN <- mapvalues(scna_data$CN, c(-2, -1, 0, 1, 2), c("DeepDel", "Del",
                                                              "wt", "ShallowAmp",
                                                              "Amp"))
scna_data <- scna_data[scna_data$CN != "wt",]
colnames(scna_data) <- c("Sample_name", "Gene", "CN")
scna_data <- scna_data[c("Gene", "Sample_name", "CN")]

#scna_data <- scna_data[scna_data$Gene %in% unique(oncoprint_data$Hugo_Symbol),]

# filter scna and snv data to only contain lipocalins
scna_data <- scna_data[tolower(scna_data$Gene) %in% tolower(lipocalin_genes),]
oncoprint_data_melanoma <- oncoprint_data_melanoma[tolower(oncoprint_data_melanoma$Hugo_Symbol) %in% tolower(lipocalin_genes),]

# filter Del and shallowAmp from scna data
scna_data_filt <- scna_data[!(scna_data$CN %in% c("Del", "ShallowAmp")),]

maf_object <- read.maf(maf = oncoprint_data, cnTable = scna_data)
maf_object_filt <- read.maf(maf = oncoprint_data, cnTable = scna_data_filt)

pdf("oncoprint_lipocalins_full_melanoma.pdf")
oncoplot(maf_object, includeColBarCN = TRUE, colors = vc_cols, top=50)
dev.off()

pdf("oncoprint_lipocalins_deep_events_melanoma.pdf")
oncoplot(maf_object_filt, includeColBarCN = TRUE, colors = vc_cols, top=50)
dev.off()

scna_data <- read.csv("../../../datasets/cnv_data/prostate_cna_obp_data.csv", header = T)
scna_data$X <- str_replace_all(scna_data$X, "\\.", "-")

oncoprint_data_prostate <- oncoprint_data[oncoprint_data$Tumor_Sample_Barcode %in% scna_data$X,] 

scna_data <- scna_data %>% pivot_longer(names_to = "Gene", values_to = "CN", cols = -X)

scna_data$CN <- mapvalues(scna_data$CN, c(-2, -1, 0, 1, 2), c("DeepDel", "Del",
                                                              "wt", "ShallowAmp",
                                                              "Amp"))
scna_data <- scna_data[scna_data$CN != "wt",]
colnames(scna_data) <- c("Sample_name", "Gene", "CN")
scna_data <- scna_data[c("Gene", "Sample_name", "CN")]

#scna_data <- scna_data[scna_data$Gene %in% unique(oncoprint_data$Hugo_Symbol),]

# filter scna and snv data to only contain lipocalins
scna_data <- scna_data[tolower(scna_data$Gene) %in% tolower(lipocalin_genes),]
oncoprint_data_prostate <- oncoprint_data_prostate[tolower(oncoprint_data_prostate$Hugo_Symbol) %in% tolower(lipocalin_genes),]

# filter Del and shallowAmp from scna data
scna_data_filt <- scna_data[!(scna_data$CN %in% c("Del", "ShallowAmp")),]

maf_object <- read.maf(maf = oncoprint_data, cnTable = scna_data)
maf_object_filt <- read.maf(maf = oncoprint_data, cnTable = scna_data_filt)

pdf("oncoprint_lipocalins_full_prostate.pdf")
oncoplot(maf_object, includeColBarCN = TRUE, colors = vc_cols, top=50)
dev.off()

pdf("oncoprint_lipocalins_deep_events_prostate.pdf")
oncoplot(maf_object_filt, includeColBarCN = TRUE, colors = vc_cols, top=50)
dev.off()

scna_data <- read.csv("../../../datasets/cnv_data/uterine_cna_obp_data.csv", header = T)
scna_data$X <- str_replace_all(scna_data$X, "\\.", "-")

oncoprint_data_uterine <- oncoprint_data[oncoprint_data$Tumor_Sample_Barcode %in% scna_data$X,] 

scna_data <- scna_data %>% pivot_longer(names_to = "Gene", values_to = "CN", cols = -X)

scna_data$CN <- mapvalues(scna_data$CN, c(-2, -1, 0, 1, 2), c("DeepDel", "Del",
                                                              "wt", "ShallowAmp",
                                                              "Amp"))
scna_data <- scna_data[scna_data$CN != "wt",]
colnames(scna_data) <- c("Sample_name", "Gene", "CN")
scna_data <- scna_data[c("Gene", "Sample_name", "CN")]

#scna_data <- scna_data[scna_data$Gene %in% unique(oncoprint_data$Hugo_Symbol),]

# filter scna and snv data to only contain lipocalins
scna_data <- scna_data[tolower(scna_data$Gene) %in% tolower(lipocalin_genes),]
oncoprint_data_uterine <- oncoprint_data_uterine[tolower(oncoprint_data_uterine$Hugo_Symbol) %in% tolower(lipocalin_genes),]

# filter Del and shallowAmp from scna data
scna_data_filt <- scna_data[!(scna_data$CN %in% c("Del", "ShallowAmp")),]

maf_object <- read.maf(maf = oncoprint_data, cnTable = scna_data)
maf_object_filt <- read.maf(maf = oncoprint_data, cnTable = scna_data_filt)

pdf("oncoprint_lipocalins_full_uterine.pdf")
oncoplot(maf_object, includeColBarCN = TRUE, colors = vc_cols, top=50)
dev.off()

pdf("oncoprint_lipocalins_deep_events_uterine.pdf")
oncoplot(maf_object_filt, includeColBarCN = TRUE, colors = vc_cols, top=50)
dev.off()


