# ================================
# CARGAR LIBRERÍAS
# ================================
library(readxl)
library(dplyr)

# ================================
# CARGAR TABLE S1
# ================================
tabla <- read_excel("mmc1.xlsx")

# Ver columnas
colnames(tabla)

# ================================
# LISTA DE GENES
# ================================
genes <- c("INTS2","INTS8","MED12","MED19","MED30",
           "SUPT5H","SUPT6H","PAF1","CTR9","POLR2B",
           "HSPA9","PHB","PHB2","DNAJC19","TIMM23B",
           "NCBP2","SRRT","DDX41","GATA1","GAB2")

# ================================
# FILTRAR GENES
# ================================
tabla_filtrada <- tabla %>%
  filter(gene %in% genes)

# ================================
# EXTRAER GUÍAS
# ================================
guias_exp <- tabla_filtrada %>%
  select(
    gene,
    transcript,
    sgID_A,
    sgID_B,
    `targeting sequence A`,
    `targeting sequence B`
  )

# ================================
# LIMPIAR NOMBRES
# ================================
colnames(guias_exp) <- c(
  "gene",
  "transcript",
  "sgRNA_A_ID",
  "sgRNA_B_ID",
  "sgRNA_A_seq",
  "sgRNA_B_seq"
)

# ================================
# GUARDAR CSV
# ================================
write.csv(guias_exp, "guias_experimentales_referencia.csv", row.names = FALSE)

# ================================
# CHECK RÁPIDO
# ================================
print(head(guias_exp))
print(unique(guias_exp$gene))
print(n_distinct(guias_exp$gene))
