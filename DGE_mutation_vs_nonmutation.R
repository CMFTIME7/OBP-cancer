library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(dplyr)
library(limma)
library(edgeR)

# Function to set parameters based on cancer type
setParameters <- function(cancer_type) {
  switch(cancer_type,
         "Ovarian cancer" = list(
           paraStudy = "GTEX",
           paraPrimarySiteGTEx = "Ovary",
           paraPrimaryTissueGTEx = "Ovary",
           paraCohort = "TCGA Ovarian Cancer",
           paraPrimarySiteTCGA = "Ovary",
           paraDatasets = "TCGA.OV.sampleMap/OV_clinicalMatrix",
           markerGenes = c("CLDN4", "HLA-DRA", "CLDN3"),
           excludeSamples = "TCGA-25-1870-01",
           paraxlim = c(-5, 13),
           paraylim = c(0, 150)
         ),
         "Breast cancer" = list(
           paraStudy = "GTEX",
           paraPrimarySiteGTEx = "Breast",
           paraPrimaryTissueGTEx = "Breast - Mammary Tissue",
           paraCohort = "TCGA Breast Cancer",
           paraPrimarySiteTCGA = "Breast",
           paraDatasets = "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
           markerGenes = c("ERBB2", "KRT19", "KRT18"),
           excludeSamples = NA,
           paraxlim = c(-5, 5),
           paraylim = c(0, 150)
         ),
         "Lung cancer" = list(
           paraStudy = "GTEX",
           paraPrimarySiteGTEx = "Lung",
           paraPrimaryTissueGTEx = "Lung",
           paraCohort = "TCGA Lung Cancer",
           paraPrimarySiteTCGA = "Lung",
           paraDatasets = "TCGA.LUNG.sampleMap/LUNG_clinicalMatrix",
           markerGenes = c("SPAG5", "KIF23", "RAD54L"),
           excludeSamples = "GTEX-SUCS-0626-SM-5CHQE",
           paraxlim = c(-4, 5),
           paraylim = c(0, 75)
         ),
         "Uterine cancer" = list(
           paraStudy = "GTEX",
           paraPrimarySiteGTEx = "Uterus",
           paraPrimaryTissueGTEx = "Uterus",
           paraCohort = "TCGA Uterine Carcinosarcoma",
           paraPrimarySiteTCGA = "Uterus",
           paraDatasets = "TCGA.UCS.sampleMap/UCS_clinicalMatrix",
           markerGenes = c("ESPL1", "PTTG1", "PRAME"),
           excludeSamples = NA,
           paraxlim = c(-5, 10),
           paraylim = c(0, 75)
         ),
         "Colorectal cancer" = list(
           paraStudy = "GTEX",
           paraPrimarySiteGTEx = "Colon",
           paraPrimaryTissueGTEx = "Colon",
           paraCohort = "TCGA Colon Cancer",
           paraPrimarySiteTCGA = "Colon",
           paraDatasets = "TCGA.COAD.sampleMap/COAD_clinicalMatrix",
           markerGenes = c("KLK6", "SFTA2", "LEMD1"),
           excludeSamples = c("GTEX-SUCS-1026-SM-5CHTC",
                              "GTEX-WFG8-1526-SM-5CHSI",
                              "GTEX-WFON-1426-SM-5CHT1",
                              "GTEX-WFG8-1726-SM-5CHTE"),
           paraxlim = c(-8, 8),
           paraylim = c(0, 120)
         ),
         "Prostate cancer" = list(
           paraStudy = "GTEX",
           paraPrimarySiteGTEx = "Prostate",
           paraPrimaryTissueGTEx = "Prostate",
           paraCohort = "TCGA Prostate Cancer",
           paraPrimarySiteTCGA = "Prostate",
           paraDatasets = "TCGA.PRAD.sampleMap/PRAD_clinicalMatrix",
           markerGenes = c("EPHA10", "HPN", "SLC4A2"),
           excludeSamples = NA,
           paraxlim = c(-8, 8),
           paraylim = c(0, 100)
         ),
         "Melanoma" = list(
           paraStudy = "GTEX",
           paraPrimarySiteGTEx = "Skin",
           paraPrimaryTissueGTEx = "Skin",
           paraCohort = "TCGA Melanoma",
           paraPrimarySiteTCGA = "Skin",
           paraDatasets = "TCGA.SKCM.sampleMap/SKCM_clinicalMatrix",
           markerGenes = NA,
           excludeSamples = NA,
           paraxlim = c(-8, 8),
           paraylim = c(0, 100)
         )
  )
}


# Set parameters for the specific cancer type
cancer_type = "Ovarian cancer"
params <- setParameters(cancer_type)


# Step 4
data(XenaData)
write.csv(XenaData, "00_tblXneaHubInfo.csv")

# Step 5-a
GeneExpectedCnt_toil = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count")
XenaQuery(GeneExpectedCnt_toil) %>%
  XenaDownload(destdir = "./")

# Step 5-b
# paraCohort = "TCGA Breast Cancer"
# paraDatasets = "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"

Clin_TCGA = XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterCohorts = params$paraCohort) %>%
  XenaFilter(filterDatasets = params$paraDatasets)
XenaQuery(Clin_TCGA) %>%
  XenaDownload(destdir = "./")

# Step 5-c
Surv_TCGA = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TCGA_survival_data")
XenaQuery(Surv_TCGA) %>%
  XenaDownload(destdir = "./")

# Step 5-d
Pheno_GTEx = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGTEx_phenotype")
XenaQuery(Pheno_GTEx) %>%
  XenaDownload(destdir = "./")

# Step 6-a
filterGTEx01 = fread("TcgaTargetGtex_phenotype.txt.gz")
names(filterGTEx01) = gsub("\\_", "", names(filterGTEx01))

# paraStudy = "GTEX"
# paraPrimarySiteGTEx = "Breast"
# paraPrimaryTissueGTEx = "Breast - Mammary Tissue"

filterGTEx02 = subset(filterGTEx01,
                      study == params$paraStudy &
                        primarysite == params$paraPrimarySiteGTEx &
                        grepl(params$paraPrimaryTissueGTEx, filterGTEx01$"primary disease or tissue"))

# Step 6-b
filterTCGA01 = fread(params$paraDatasets)
names(filterTCGA01) = gsub("\\_", "", names(filterTCGA01))

paraSampleType = "Primary Tumor"
# paraPrimarySiteTCGA = "Breast"



filterTCGA02 = subset(filterTCGA01,
                      sampletype == paraSampleType &
                        primarysite == params$paraPrimarySiteTCGA)

# Select sampleIDs with missense mutations
mutation_samples = read.csv("Missense_mutations_sampleIDs.csv")

# Extract a single column
col_names <- params$paraPrimarySiteTCGA
mutation_sampleIDs <- mutation_samples[[params$paraPrimarySiteTCGA]]

filterTCGA02_m = subset(filterTCGA02,
                      sampleID %in% mutation_sampleIDs)

# Select sampleIDs without mutations

filterTCGA_nm = subset(filterTCGA02, !(sampleID %in% mutation_sampleIDs))


# Step 6-c
filterExpr = c(filterTCGA_nm$sampleID, filterTCGA02_m$sampleID, "sample")
ExprSubsetBySamp = fread("TcgaTargetGtex_gene_expected_count.gz", select = filterExpr)

# Check for nm and m sampleID respectively

filterExpr_nm = c(filterTCGA_nm$sampleID, "sample")
filterExpr_m = c(filterTCGA02_m$sampleID, "sample")

ExprSubsetBySamp_nm = fread("TcgaTargetGtex_gene_expected_count.gz", select = filterExpr_nm)
ExprSubsetBySamp_m = fread("TcgaTargetGtex_gene_expected_count.gz", select = filterExpr_m)

# Change the first 4 letters

# Replace 'TCGA' with 'TCGB' in the column names
colnames(ExprSubsetBySamp_nm) <- gsub("^TCGA", "TCGB", colnames(ExprSubsetBySamp_nm))
ExprSubsetBySamp = merge(ExprSubsetBySamp_nm, ExprSubsetBySamp_m, by = "sample") 
# cbind(ExprSubsetBySamp_nm, ExprSubsetBySamp_m)

# Step 7
probemap = fread("zz_gencode.v23.annotation.csv", select = c(1, 2))
exprALL = merge(probemap, ExprSubsetBySamp, by.x = "id", by.y = "sample")
genesPC = fread("zz_gene.protein.coding.csv")
exprPC = subset(exprALL, gene %in% genesPC$Gene_Symbol)

exprFinal = exprPC[!(duplicated(exprPC$gene) |
                       duplicated(exprPC$gene, fromLast = T))] # Check this for genes

# Step 8
write.csv(exprFinal, "00_ExpectedCnt.csv")

# Step 12
exprBT = exprFinal[, -c(1:2)] #왜 3까지였지?
rownames(exprBT) = exprFinal$gene
exprBT = round(((2^exprBT)-1), 0)
write.csv(exprBT, "01_EpectedCntBT.csv")

# Step 13
# exprLIMMA = exprBT[, -params$excludeSamples]

if (!is.na(params$excludeSamples) && length(params$excludeSamples) > 0) {
  # Check if excludeSamples are in the column names
  if (all(params$excludeSamples %in% colnames(exprBT))) {
    # Create a logical vector to keep columns not in excludeSamples
    samplesToKeep = !colnames(exprBT) %in% params$excludeSamples
    exprLIMMA = exprBT[, .SD, .SDcols = samplesToKeep]
  } else {
    # If some of the excludeSamples are not in the column names, handle this case
    warning("Some excludeSamples not found in column names.")
    exprLIMMA = exprBT
  }
} else {
  # If excludeSamples is NA or empty, use all samples
  exprLIMMA = exprBT
}

x = DGEList(exprLIMMA)
# x = DGEList(exprBT)

# Step 14
snames = colnames(x)
group = substr(snames, 1, 4)
x$samples$group = group

# Step 16
keep.exprs = filterByExpr(x, group = group, min.count = 1)
x = x[keep.exprs, , keep.lib.sizes = FALSE]

# Step 17
x  = calcNormFactors(x, method = "upperquartile")

# Step 18
design = model.matrix(~0 + group)

# Step 19
colnames(design) = gsub("group", "", colnames(design))
contr.matrix = makeContrasts(TCGAvsTCGB = TCGA - TCGB,
                             levels = colnames(design))

# Step 20
v = voom(x, design)

# Step 21
vfit = lmFit(v, design)
vfit = contrasts.fit(vfit, contrasts = contr.matrix)

# Step 22
efit = eBayes(vfit)

# Step 23
tfit = treat(vfit, lfc = 0.58)

# Step 24
DGEsTreat = topTreat(tfit, n = Inf)
geneIDs = as.numeric(rownames(DGEsTreat))
genesDGE = exprFinal[geneIDs, "gene"]
rownames(DGEsTreat) = genesDGE$gene

write.csv(DGEsTreat, "01_DGEsTreat.csv")

# Step 25
voomExpr = v$E
write.csv(voomExpr, "01_voomExpr.csv")

# Script for generating a volcano plot
library(EnhancedVolcano)
library(gtable)
library(gridExtra)
library(ggalt)

dataExpr = read.csv("01_DGEsTreat.csv")

# lipocalins  = c("OBP2A", "OBP2B", "LCN2", "LCN1", "PAEP", "AMBP", "APOD",
#                 "APOM", "CRABP2", "RBP4", "LCN9", "LCN10", "LCN15")

lipocalins = c("LCN1","LCN2","LCN6","LCN8","LCN9","LCN10","LCN12","LCN15","OBP2A","OBP2B",
               "AMBP","APOD","APOM","C8G","ORM1","ORM2","PAEP","PTGDS","RBP4")

keyvals = ifelse(dataExpr$X %in% lipocalins, 
                 "black", 
                 ifelse(dataExpr$X %in% params$markerGenes,
                        "red4",
                        ifelse(dataExpr$logFC < -0.58 & 
                                 dataExpr$adj.P.Val < 0.05, 
                               "blue", 
                               ifelse(dataExpr$logFC > 0.58 & 
                                        dataExpr$adj.P.Val < 0.05, 
                                      "red", 
                                      "gray")))
)

keyvals[is.na(keyvals)] = "gray"
names(keyvals)[keyvals == "red"] = "Up-regulated"
names(keyvals)[keyvals == "blue"] = "Down-regulated"
names(keyvals)[keyvals == "gray"] = "Not significant"
names(keyvals)[keyvals == "red4"] = "Cancer markers"
names(keyvals)[keyvals == "black"] = "Lipocalins"


keyvals.shape = ifelse(dataExpr$X == "OBP2A", 18, 
                       ifelse(dataExpr$X == "OBP2B", 17, 
                              ifelse(dataExpr$X == "LCN2", 15, 
                                     ifelse(dataExpr$X %in% params$markerGenes,
                                            16,
                                            ifelse(
                                              dataExpr$X %in% lipocalins,
                                              42,
                                              1)))))

keyvals.size = ifelse(dataExpr$X == "OBP2A", 3, 
                      ifelse(dataExpr$X == "OBP2B", 3, 
                             ifelse(dataExpr$X == "LCN2", 3, 
                                    ifelse(dataExpr$X %in% params$markerGenes, 3,
                                           ifelse(dataExpr$X %in% lipocalins, 6, 1)))))


names(keyvals.shape)[keyvals.shape == 18] = "OBP2A"
names(keyvals.shape)[keyvals.shape == 17] = "OBP2B"
names(keyvals.shape)[keyvals.shape == 15] = "LCN2"
names(keyvals.shape)[keyvals.shape == 16] = "Cancer markers"
names(keyvals.shape)[keyvals.shape == 42] = "Other lipocalins"
names(keyvals.shape)[keyvals.shape == 1] = "Other proteins"

p = EnhancedVolcano(
  dataExpr,
  lab = dataExpr$X,
  x = "logFC",
  y = "adj.P.Val",
  selectLab = lipocalins,
  xlim = params$paraxlim,
  ylim = c(0,50),
  xlab = "log2 fold change",
  ylab = "-log(FDR)",
  axisLabSize = 12,
  pCutoff = 0.05,
  FCcutoff = 0.58,
  # title = paste0(paste0("TCGA Tumor vs GTEx Normal (", cancer_type),")"),
  title = paste0(paste0("TCGA Tumor mutation vs non-mutation (", cancer_type),")"),
  titleLabSize = 14,
  subtitle = NULL,
  caption = NULL,
  pointSize = keyvals.size,
  labSize = 3,
  labCol = 'black',
  #labFace = 'bold',
  boxedLabels = T,
  colAlpha = 4/5,
  colCustom = keyvals,
  shapeCustom = keyvals.shape,
  legendPosition = "right", 
  legendLabSize = 8,
  legendIconSize = 3,
  drawConnectors = TRUE,
  arrowheads = FALSE,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
)

p