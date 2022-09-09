library(magrittr)
library(matrixStats)
library(Matrix)
library(do)

# (Mus Spermatogenesis) GSE107644 -----
path <- "../../data/(Mus Spermatogenesis) GSE107644/GEO Raw Data/GSE107644_RAW/"
filename <- list.files(path)
testis_scs <- as.data.frame(array(dim = c(23154,0)))
for (i in 1:length(filename)) {
  tmp_dir <- paste0(path,filename[i])
  tmp_read <- read.table(tmp_dir,header = TRUE, stringsAsFactors = FALSE)
  testis_scs <- data.frame(testis_scs,tmp_read[,-1])
}
testis_scs <- cbind(GeneName = tmp_read[,1],testis_scs)

cell_anno <- read.table("../../data/(Mus Spermatogenesis) GSE107644/GEO Raw Data/GSE107644_barcode_information.txt",header = TRUE)
celltype <- c()
for (i in 1:nrow(cell_anno)) {
  tmp_celltype <- data.frame(unlist(strsplit(cell_anno[i,3],split = "_")))[2,1]
  celltype <- c(celltype,tmp_celltype)
}
celltype <- gsub("TypeBS","BS",celltype)
celltype <- gsub("TypeBG2M","BG2",celltype)
celltype <- gsub("RS1o2","RS2",celltype)
celltype <- gsub("RS3o4","RS4",celltype)
celltype <- gsub("RS5o6","RS6",celltype)
celltype <- gsub("RS7o8","RS8",celltype)
cell_anno$celltype <- celltype

celltype_list <- c("A1","ln","BS","BG2","G1","ePL","mPL","lPL","L","Z","eP","mP",
  "lP","D","MI","MII","RS2","RS4","RS6","RS8")
testis_scs_median <- data.frame(GeneName = testis_scs$GeneName)
for (celltype in celltype_list) {
  cell.idx <- cell_anno[cell_anno$celltype == celltype,3]
  cell.data <- testis_scs[,names(testis_scs) %in% cell.idx] %>% as.matrix()
  testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
  names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Mouse",names(testis_scs_median)[2:ncol(testis_scs_median)],sep = "_")

gse107644 <- testis_scs_median
names(gse107644) <- c("GeneName",paste("GSE107644",names(gse107644)[2:ncol(gse107644)],sep = "_"))

# (Hom Spermatogenesis) GSE106487 -----
path <- "../../data/(Hom Spermatogenesis) GSE106487/GEO Raw Data/GSE106487_RAW/"
filename <- list.files(path)
testis_scs <- as.data.frame(array(dim = c(24153,0)))
for (i in 1:length(filename)) {
  tmp_dir <- paste0(path,filename[i])
  tmp_read <- read.table(tmp_dir,header = TRUE, stringsAsFactors = FALSE)
  testis_scs <- data.frame(testis_scs,tmp_read[,-1])
}
testis_scs <- cbind(Symbol = tmp_read[,1],testis_scs)
names(testis_scs)[1] <- "GeneName"
names(testis_scs) <- gsub("X","",names(testis_scs))

celltype_anno <- read.csv("../../data/(Hom Spermatogenesis) GSE106487/GEO Raw Data/annotation_celltype.csv")
cell_list <- data.frame(cell = names(testis_scs)[-1])
cell_list$cell <- gsub("X","",cell_list$cell)
cell_anno <- merge(cell_list,celltype_anno,all.x = TRUE) %>% 
  .[!is.na(.$class),]

celltype_list <- c("SSC.subpopulation1","SSC.subpopulation2","Diff.ing.SPG",
  "Diff.ed.SPG","L1","L2","L3","Z","P","D","SPC7","S1","S2","S3","S4","ST")
testis_scs_median <- data.frame(GeneName = testis_scs$GeneName)
for (celltype in celltype_list) {
  cell.idx <- cell_anno[cell_anno$class == celltype,1]
  cell.data <- testis_scs[,names(testis_scs) %in% cell.idx] %>% as.matrix()
  testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
  names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Human",names(testis_scs_median)[2:ncol(testis_scs_median)],sep = "_")

gse106487 <- testis_scs_median
names(gse106487) <- c("GeneName",paste("GSE106487",names(gse106487)[2:ncol(gse106487)],sep = "_"))

# (Hom Somat) GSE124263 -----
path <- "../../data/(Hom Somat) GSE124263/GEO Raw Data/GSE124263_RAW/"
prefix_list <- c(
  "GSM3526583_d2_I_",
  "GSM3526584_D2Total_",
  "GSM3526585_d7_I_",
  "GSM3526586_D7total_",
  "GSM3526587_A1_I_",
  "GSM3526588_A1Total_", 
  "GSM3526589_A2_I_",
  "GSM3526590_A2_total_")
testis_scs <- as.data.frame(array(dim = c(32738,0)))
for (prefix in prefix_list) {
  temp_path <- paste0(path,prefix)
  
  barcode.path <- paste0(temp_path,"barcodes.tsv")
  genes.path <- paste0(temp_path,"genes.tsv")
  matrix.path <- paste0(temp_path,"matrix.mtx")
  
  data <- readMM(file = matrix.path)
  gene.names <- read.delim(file = genes.path,header = FALSE)
  barcode.names <- read.delim(barcode.path,header = FALSE)
  
  rownames(data) = gene.names$V1
  colnames(data) = barcode.names$V1
  
  data <- as.matrix(data)
  
  testis_scs <- cbind(testis_scs,data)
}
testis_scs_names <- gsub("-1","",names(testis_scs))
names(testis_scs) <- testis_scs_names

cell_anno <- read.csv("../../data/(Hom Somat) GSE124263/GEO Raw Data/annotation_celltype.csv")
cell_anno$class <- Replace(
  cell_anno$class,
  pattern = c(
    "[(].*[)]:",
    " $:",
    "-\\d:",
    "\\d:",
    "Leydig Cells:LC",
    "Peri-Tubular Myoid:PTM",
    "EC/Blood:EC",
    "Spermatogonia:SPG",
    "Spermatid:St",
    "Spermatocytes:SPC",
    "Macrophage:M",
    "Transition Cell:TC"))
sample <- c()
barcode <- c()
for (i in 1:nrow(cell_anno)) {
  temp_strsp <- as.data.frame(unlist(strsplit(cell_anno[i,1],split = "_")))
  sample <- c(sample,temp_strsp[1,])
  barcode <- c(barcode,temp_strsp[2,])
}
cell_anno$sample <- sample
cell_anno$barcode <- barcode

celltype_list <- c("SC","LC","PTM","GC","EC","PreSPG","PGCL","SPG","St","SPC","M",
  "SSC","TC","Diff SPG","Early diff SPG")
prenatal_sample_list <- c("D2T","D2I","D7T","D7I")
posnatal_sample_list <- c("A1T","A1I","A2T","A2I")
testis_scs_median <- data.frame(Ensembl.ID = rownames(testis_scs))
for (celltype in celltype_list) {
  prenatal.cell.idx <- cell_anno[(cell_anno$class == celltype) & (cell_anno$sample %in% prenatal_sample_list),4]
  posnatal.cell.idx <- cell_anno[(cell_anno$class == celltype) & (cell_anno$sample %in% posnatal_sample_list),4]
  
  prenatal.cell.data <- testis_scs[,names(testis_scs) %in% prenatal.cell.idx] %>% as.matrix()
  posnatal.cell.data <- testis_scs[,names(testis_scs) %in% posnatal.cell.idx] %>% as.matrix()
  
  if (!ncol(prenatal.cell.data) == 0) {
    testis_scs_median <- cbind(testis_scs_median,rowMedians(prenatal.cell.data))
    names(testis_scs_median)[ncol(testis_scs_median)] <- paste("Prenatal",celltype,sep = "_")
  }
  if (!ncol(posnatal.cell.data) == 0) {
    testis_scs_median <- cbind(testis_scs_median,rowMedians(posnatal.cell.data))
    names(testis_scs_median)[ncol(testis_scs_median)] <- paste("Posnatal",celltype,sep = "_")
  }
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Human",names(testis_scs_median)[2:ncol(testis_scs_median)],sep = "_")

gse124263 <- testis_scs_median
names(gse124263) <- c("Ensembl.ID",paste("GSE124263",names(gse124263)[2:ncol(gse124263)],sep = "_"))

# (Hom Gonad) GSE86146 -----
path <- "../../data/(Hom Gonad) GSE86146/GEO Raw Data/GSE86146_RAW/"
filename <- list.files(path)
testis_scs <- as.data.frame(array(dim = c(24153,0)))
for (i in 1:length(filename)) {
  tmp_dir <- paste0(path,filename[i])
  tmp_read <- read.table(tmp_dir,header = TRUE, stringsAsFactors = FALSE)
  testis_scs <- data.frame(testis_scs,tmp_read[,-1])
}
testis_scs <- cbind(GeneName = tmp_read[,1],testis_scs)

cell_anno <- read.csv("../../data/(Hom Gonad) GSE86146/GEO Raw Data/annotation_celltype.csv")
sample <- c()
for (i in 1:nrow(cell_anno)) {
  temp_strsp <- as.data.frame(unlist(strsplit(cell_anno[i,1],split = "_")))
  sample <- c(sample,temp_strsp[2,])
}
cell_anno$sample <- sample

celltype_list <- c(
  "Female_FGC#1","Female_FGC#2","Female_FGC#3","Female_FGC#4",
  "Female_Soma#1","Female_Soma#2","Female_Soma#3","Female_Soma#4",
  "Male_FGC#1","Male_FGC#2","Male_FGC#3","Male_Soma#1","Male_Soma#2","Male_Soma#3","Male_Soma#4")
sample_list <- c("4W","5W","7W","8W","10W","11W","12W","14W","18W","20W","23W","24W","26W","25W","9W","19W","21W")
testis_scs_median <- data.frame(GeneName = testis_scs$GeneName)
for (celltype in celltype_list) {
  for (time in sample_list) {
    cell.idx <- cell_anno[(cell_anno$class == celltype) & (cell_anno$sample == time),1]
    cell.data <- testis_scs[,names(testis_scs) %in% cell.idx] %>% as.matrix()
    if (!ncol(cell.data) == 0) {
      testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
      names(testis_scs_median)[ncol(testis_scs_median)] <- paste(time,celltype,sep = "_")
    }
  }
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Human",names(testis_scs_median)[2:ncol(testis_scs_median)],sep = "_")

gse86146 <- testis_scs_median
names(gse86146) <- c("GeneName",paste("GSE86146",names(gse86146)[2:ncol(gse86146)],sep = "_"))

# (Hom Adult SPG) GSE108977 -----
testis_scs <- read.table("../../data/(Hom Adult SPG) GSE108977/GSE108977_HumanSpg-C1_Gematrix.txt",header = TRUE)
celltype_list <- paste("HumanSpg.Seq",c(1:7,9:10),sep = "")
testis_scs_median <- data.frame(GeneName = testis_scs$ID)
for (celltype in celltype_list) {
  cell.idx <- grep(paste(celltype,"_",sep = ""),names(testis_scs))
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
  names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Human_Adult_Spg.Seq",c(1:7,9:10),sep = "")

gse108977 <- testis_scs_median
names(gse108977) <- c("GeneName",paste("GSE108977",names(gse108977)[2:ncol(gse108977)],sep = "_"))

# (Mus Adult SPG) GSE108974 -----
testis_scs <- read.table("../../data/(Mus Adult SPG) GSE108974/GSE108974_Ad-Id4GFP-C1_Gematrix.txt",header = TRUE)
celltype_list <- paste("AdId4eGFPSeq",1:4,sep = "")
testis_scs_median <- data.frame(GeneName = testis_scs$ID)
for (celltype in celltype_list) {
  cell.idx <- grep(paste(celltype,"_",sep = ""),names(testis_scs))
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
  names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Mouse_Adult_Spg.Seq",1:4,sep = "")

gse108974 <- testis_scs_median
names(gse108974) <- c("GeneName",paste("GSE108974",names(gse108974)[2:ncol(gse108974)],sep = "_"))

# (Mus P6 SPG) GSE108970 -----
testis_scs <- read.table("../../data/(Mus P6 SPG) GSE108970/GSE108970_P6-Id4GFP-C1_Gematrix.txt",header = TRUE)
celltype_list <- paste("Id4GFP.Seq",1:5,sep = "")
testis_scs_median <- data.frame(GeneName = testis_scs$ID)
for (celltype in celltype_list) {
  cell.idx <- grep(paste(celltype,"_",sep = ""),names(testis_scs))
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
  names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Mouse_P6_Spg.Seq",1:5,sep = "")

gse108970 <- testis_scs_median
names(gse108970) <- c("GeneName",paste("GSE108970",names(gse108970)[2:ncol(gse108970)],sep = "_"))

# (Mus P3&7 SPG) GSE82174 -----
testis_scs <- read.table("../../data/(Mus P3&7 SPG) GSE82174/GSE82174_TPM_Matrix.txt",header = TRUE)
celltype_list <- c("WTp7A","WTp7B","WTp3A","WTp3B")
testis_scs_median <- data.frame(GeneName = testis_scs$Gene)
for (celltype in celltype_list) {
  cell.idx <- grep(paste(celltype,"_",sep = ""),names(testis_scs))
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
  names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- c("Mouse_P7_Spg.A","Mouse_P7_Spg.B","Mouse_P3_Spg.A","Mouse_P3_Spg.B")

gse82174 <- testis_scs_median
names(gse82174) <- c("GeneName",paste("GSE82174",names(gse82174)[2:ncol(gse82174)],sep = "_"))

# (Mus E10.5-E16.5 Gonad) GSE97519 -----
testis_scs <- read.table("../../data/(Mus E10.5-E16.5 Gonad) GSE97519/GSE97519_XY_NR5A1-eGFP_single-cell_fetal_gonads.txt",header = TRUE)
celltype_list <- c(
  "E10.5_XY_20140428",
  "E11.5_XY_20141103","E11.5_XY_20150107",
  "E12.5_XY_20140526","E12.5_XY_20141210",
  "E13.5_XY_20130918","E13.5_XY_20140528",
  "E16.5_XY_20150202","E16.5_XY_20150223")
testis_scs_median <- data.frame(GeneName = rownames(testis_scs))
for (celltype in celltype_list) {
  cell.idx <- grep(paste(celltype,"_",sep = ""),names(testis_scs))
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
  names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- c(
  "Mouse_E10.5_GN.Seq1",
  "Mouse_E11.5_GN.Seq1","Mouse_E11.5_GN.Seq2",
  "Mouse_E12.5_GN.Seq1","Mouse_E12.5_GN.Seq2",
  "Mouse_E13.5_GN.Seq1","Mouse_E13.5_GN.Seq2",
  "Mouse_E16.5_GN.Seq1","Mouse_E16.5_GN.Seq2")

gse97519 <- testis_scs_median
names(gse97519) <- c("GeneName",paste("GSE97519",names(gse97519)[2:ncol(gse97519)],sep = "_"))

# (Mus P5.5 SCS) GSE107711 -----
testis_scs <- read.table("../../data/(Mus P5.5 SCS) GSE107711/GSE107711_Salmon.Expected.count.matrix.txt",sep = "\t",skip = 2,header = TRUE)
celltype_list <- paste("Cluster",c("I","II","III","IV"),sep = ".")
testis_scs_median <- data.frame(GeneName = testis_scs$Group)
for (celltype in celltype_list) {
  cell.idx <- grep(paste(celltype,".",sep = ""),names(testis_scs))
  cell.data <- testis_scs[,cell.idx ] %>% as.matrix()
  testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
  names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- c("Mouse_P5.5_Gonocyte","Mouse_P5.5_PrimitiveSC","Mouse_P5.5_Progenitor","Mouse_P5.5_DiffedPrimeS")

gse107711 <- testis_scs_median
names(gse107711) <- c("GeneName",paste("GSE107711",names(gse107711)[2:ncol(gse107711)],sep = "_"))

# (Hom Somat) GSE112013 -----
testis_scs <- read.table("../../data/(Hom Somat) GSE112013/GSE112013_Combined_UMI_table.txt",header = TRUE)
celltype_anno <- read.table("../../data/(Hom Somat) GSE112013/41422_2018_99_MOESM9_ESM.txt",header = TRUE)
celltype_list <- 1:13
testis_scs_median <- data.frame(GeneName = testis_scs$Gene)
for (celltype in celltype_list) {
  cell.idx <- which(names(testis_scs) %in% gsub("-",".",celltype_anno[celltype_anno$Final_clusters == celltype,1]))
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
  names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste(
  "Human",
  c("SSCs","DifferentiatingSpg","EarlyPrimarySPC","LatePrimarySPC","RS","ES","Sperm1","Sperm2","Macrophage","EndothelialCells","MyoidCells","SertoliCells","LeydigCells"),
  sep = "_")
  
gse112013 <- testis_scs_median
names(gse112013) <- c("GeneName",paste("GSE112013",names(gse112013)[2:ncol(gse112013)],sep = "_"))

# (Mus Somat) GSE112393 -----
testis_scs_p1 <- read.table("../../data/(Mus Somat) GSE112393/GSE112393_MergedAdultMouseST25_DGE.txt",header = TRUE,nrows = 13000,)
testis_scs_p1 <- cbind(GeneName = rownames(testis_scs_p1),testis_scs_p1)
testis_scs_p2 <- read.table("../../data/(Mus Somat) GSE112393/GSE112393_MergedAdultMouseST25_DGE.txt",header = FALSE,skip = 13001)
names(testis_scs_p2) <- names(testis_scs_p1)

testis_scs <- rbind(testis_scs_p1,testis_scs_p2)
rownames(testis_scs) <- testis_scs$GeneName

celltype_anno <- read.table("../../data/(Mus Somat) GSE112393/GSE112393_MergedAdultMouseST25_PerCellAttributes.txt",header = TRUE)
celltype_list <- c(1:6,8:11)
testis_scs_median <- data.frame(GeneName = testis_scs$GeneName)
for (celltype in celltype_list) {
  cell.idx <- which(names(testis_scs) %in% celltype_anno[celltype_anno$CellType == celltype,1])
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
  names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste(
  "Mouse",
  c("InnateLymph","Macrophage","EndothelialCells","MyoidCells","LeydigCells","SertoliCells","SPG","Scytes","STids","ES"),
  sep = "_")

testis_scs_median_gc <- read.csv("../../data/(Mus Somat) GSE112393/GSE112393_MergedAdultMouseST25_12GermCellClusters_AllGeneExp.csv")
names(testis_scs_median_gc)[2:ncol(testis_scs_median_gc)] <- paste("Mouse",names(testis_scs_median_gc)[-1],sep = "_")
testis_scs_median <- merge(testis_scs_median,testis_scs_median_gc,all.x = TRUE)

gse112393 <- testis_scs_median
names(gse112393) <- c("GeneName",paste("GSE112393",names(gse112393)[2:ncol(gse112393)],sep = "_"))

# (Hom Puberty) GSE134144 -----
testis_scs <- read.table("../../data/(Hom Puberty) GSE134144/GSE134144_Pubertal_combined_UMI.txt",header = TRUE)
celltype_anno <- read.csv("../../data/(Hom Puberty) GSE134144/mmc2.csv")
celltype_anno$CellID <- paste(
  sapply(strsplit(celltype_anno$CellID,split = "-"),"[",1),
  paste(
    sapply(strsplit(celltype_anno$CellID,split = "-"),"[",3),
    sapply(strsplit(celltype_anno$CellID,split = "-"),"[",4),sep = "."),sep = ".")
celltype_anno$Age_CellType <- paste(celltype_anno$Age,celltype_anno$CellType,sep = "_")
celltype_list <- unique(celltype_anno$Age_CellType)
testis_scs_median <- data.frame(GeneName = testis_scs$gene)
for (celltype in celltype_list) {
  cell.idx <- which(names(testis_scs) %in% celltype_anno[celltype_anno$Age_CellType == celltype,1])
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  if (!length(cell.idx) == 0) {
    testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
    names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
  }
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Human",names(testis_scs_median)[2:ncol(testis_scs_median)],sep = "_")

gse134144 <- testis_scs_median
names(gse134144) <- c("GeneName",paste("GSE134144",names(gse134144)[2:ncol(gse134144)],sep = "_"))

# (Mus Gonad) GSE136220 -----
testis_scs <- read.csv("../../data/(Mus Gonad) GSE136220/GSE136220_raw_counts_matrix_pgcs_gonad_adrenal.csv")
row.names(testis_scs) <- testis_scs$index
testis_scs <- testis_scs[,-1] %>% t() %>% as.data.frame()

celltype_anno <- read.csv("../../data/(Mus Gonad) GSE136220/GSE136220_series_matrix.csv")
celltype_anno$X.Sample_title <- gsub("pgcs_","",celltype_anno$X.Sample_title)
celltype_anno$X.Sample_characteristics_ch1 <- gsub("developmental stage: ","",celltype_anno$X.Sample_characteristics_ch1)
celltype_anno$X.Sample_characteristics_ch1.1 <- gsub("Sex: ","",celltype_anno$X.Sample_characteristics_ch1.1)
celltype_anno$CellType <- paste(
  celltype_anno$X.Sample_characteristics_ch1,
  paste(
    celltype_anno$X.Sample_characteristics_ch1.1,
    celltype_anno$X.Sample_source_name_ch1,sep = "_"),sep = "_")

celltype_list <- unique(celltype_anno$CellType)
testis_scs_median <- data.frame(GeneName = rownames(testis_scs))
for (celltype in celltype_list) {
  sample.list <- celltype_anno[celltype_anno$CellType == celltype,1]
  cell.idx <- c()
  for (sample in sample.list) {
    cell.idx <- c(cell.idx,grep(sample,names(testis_scs)))
  }
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  if (!ncol(cell.data) == 0) {
    testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
    names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
  }
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Mouse",names(testis_scs_median)[2:ncol(testis_scs_median)],sep = "_")

gse136220 <- testis_scs_median
names(gse136220) <- c("Ensembl.ID",paste("GSE136220",names(gse136220)[2:ncol(gse136220)],sep = "_"))

# (Hom Somat) GSE142585 -----
testis_scs <- read.table("../../data/(Hom Somat) GSE142585/GSE142585_MergedHumanTestis4_DGE.txt",header = TRUE)
celltype_anno <- read.table("../../data/(Hom Somat) GSE142585/GSE142585_MergedHumanTestis4_PerCellAttributes.txt",header = TRUE)
celltype_anno$Barcode <- rownames(celltype_anno)

celltype_list <- unique(celltype_anno$CellType)
testis_scs_median <- data.frame(GeneName = rownames(testis_scs))
for (celltype in celltype_list) {
  cell.idx <- which(names(testis_scs) %in% celltype_anno[celltype_anno$CellType == celltype,"Barcode"])
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  if (!ncol(cell.data) == 0) {
    testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
    names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
  }
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Human",names(testis_scs_median)[2:ncol(testis_scs_median)],sep = "_")

gse142585 <- testis_scs_median
names(gse142585) <- c("GeneName",paste("GSE142585",names(gse142585)[2:ncol(gse142585)],sep = "_"))

# (Mus Gonad) GSE148032 -----
testis_scs <- read.csv("../../data/(Mus Gonad) GSE148032/GSE148032_Processeddata_UMI_counts_Male_Mouse_Germline.csv")
celltype_anno <- read.csv("../../data/(Mus Gonad) GSE148032/GSE148032_Metadata_Male_Mouse_Germline.csv")
celltype_anno <- celltype_anno[celltype_anno$genetype == "WT",]

celltype_list <- unique(celltype_anno$CellType) %>% .[-c(3,7,9)]
testis_scs_median <- data.frame(GeneName = testis_scs$X)
for (celltype in celltype_list) {
  cell.idx <- which(names(testis_scs) %in% celltype_anno[celltype_anno$CellType == celltype,3])
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  if (!ncol(cell.data) == 0) {
    testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
    names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
  }
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Mouse",names(testis_scs_median)[2:ncol(testis_scs_median)],sep = "_")

gse148032 <- testis_scs_median
names(gse148032) <- c("GeneName",paste("GSE148032",names(gse148032)[2:ncol(gse148032)],sep = "_"))

# (Hom Spermatogenesis) Lab -----
testis_scs_p1 <- readRDS("../../data/(Hom Spermatogenesis) Lab/testis_exp_MAGIC_cell1.RDS")
testis_scs_p2 <- readRDS("../../data/(Hom Spermatogenesis) Lab/testis_exp_MAGIC_cell2.RDS")
testis_scs_p3 <- readRDS("../../data/(Hom Spermatogenesis) Lab/testis_exp_MAGIC_cell3.RDS")
testis_scs_p4 <- readRDS("../../data/(Hom Spermatogenesis) Lab/testis_exp_MAGIC_cell4.RDS")

testis_scs <- Reduce(
  cbind,
  list(testis_scs_p1,testis_scs_p2,testis_scs_p3,testis_scs_p4))

celltype_anno <- read.csv("../../data/(Hom Spermatogenesis) Lab/testis_celltype.csv")

celltype_list <- unique(celltype_anno$celltype)
testis_scs_median <- data.frame(GeneName = rownames(testis_scs))
for (celltype in celltype_list) {
  cell.idx <- celltype_anno[celltype_anno$celltype == celltype,2]
  cell.data <- testis_scs[,cell.idx] %>% as.matrix()
  if (!ncol(cell.data) == 0) {
    testis_scs_median <- cbind(testis_scs_median,rowMedians(cell.data))
    names(testis_scs_median)[ncol(testis_scs_median)] <- celltype
  }
}
names(testis_scs_median)[2:ncol(testis_scs_median)] <- paste("Human",names(testis_scs_median)[2:ncol(testis_scs_median)],sep = "_")

lab.SCSeq <- testis_scs_median
names(lab.SCSeq) <- c("GeneName",paste("Lab",names(lab.SCSeq)[2:ncol(lab.SCSeq)],sep = "_"))

# Save and Load -----
load("../../data/SCSeq.Testis.RData")
save(
  gse107644,gse106487,gse124263,gse86146,gse108977,gse108974,gse108970,gse82174,
  gse97519,gse107711,gse112013,gse112393,gse134144,gse136220,gse142585,gse148032,lab.SCSeq,
  file = "../../data/SCSeq.Testis.RData")



