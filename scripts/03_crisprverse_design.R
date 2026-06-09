# ================================
# LIBRERÍAS
# ================================
library(crisprVerse)
library(crisprDesign)
library(crisprScore)
library(crisprBase)
library(AnnotationHub)
library(ensembldb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)

# ================================
# ANOTACIÓN GENÓMICA
# ================================
ah <- AnnotationHub()
edb <- ah[["AH98047"]]

data(SpCas9, package = "crisprBase")
cas9 <- SpCas9

# ================================
# LISTA DE GENES
# ================================
genes <- c("INTS2","INTS8","MED12","MED19","MED30",
           "SUPT5H","SUPT6H","PAF1","CTR9","POLR2B",
           "HSPA9","PHB","PHB2","DNAJC19","TIMM23B",
           "NCBP2","SRRT","DDX41","GATA1","GAB2")

# ================================
# BUCLE PRINCIPAL
# ================================
top3_all <- data.frame()

for (gene in genes) {
  
  cat("Procesando:", gene, "\n")
  
  # Extraer gen
  gene_gr <- genes(edb, filter = GeneNameFilter(gene))
  
  # Ajustar formato genoma
  granges <- gene_gr
  seqlevelsStyle(granges) <- "UCSC"
  genome(granges) <- "hg38"
  
  # Generar guías
  guides <- findSpacers(
    granges,
    bsgenome = BSgenome.Hsapiens.UCSC.hg38,
    crisprNuclease = cas9
  )
  
  # Calcular score
  guides <- addOnTargetScores(
    guides,
    methods = "ruleset1"
  )
  
  # Convertir a data.frame
  guides_df <- as.data.frame(guides)
  
  # Ordenar y sacar top 3
  guides_df <- guides_df[order(-guides_df$score_ruleset1), ]
  top3 <- guides_df[1:3, ]
  top3$gene <- gene
  top3_all <- rbind(top3_all, top3)
  
  # Guardar archivos por gen
  write.csv(guides_df, paste0(gene, "_crisprVerse_ALL.csv"), row.names = FALSE)
  write.csv(top3, paste0(gene, "_crisprVerse_TOP3.csv"), row.names = FALSE)
}

# ================================
# GENERAR DATASET FINAL
# ================================
write.csv(top3_all, "crisprVerse_top3_dataset.csv", row.names = FALSE)
