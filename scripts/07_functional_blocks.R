# =============================================================
# SCRIPT 07: ANÁLISIS POR BLOQUES FUNCIONALES + FIGURAS FINALES
# =============================================================
# Objetivos:
#   1. Agrupar los 20 genes en 5 bloques funcionales
#   2. Comparar Hamming y cobertura entre bloques
#   3. Generar figuras de scores SEPARADAS por herramienta
#      (corrige el error metodológico: los scores no son comparables)
#   4. Figura GC content comparativa
#
# Outputs:
#   results/tabla_bloques_funcionales.csv
#   results/tabla_gc_content.csv
#   figures/figura_scores_chopchop.png
#   figures/figura_scores_crisprverse.png
#   figures/figura_gc_content.png
#   figures/figura_bloques_hamming.png
#   figures/figura_bloques_cobertura.png
# =============================================================

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

BASE <- "/Users/jaimehache/Desktop/GENÉTICA/(iii) MÁSTER BIOINFORMÁTICA/(ii) SEGUNDO SEMESTRE/TFM/GITHUB"
setwd(BASE)
cat("📂 Directorio:", BASE, "\n\n")

# =============================================================
# 1. DEFINIR BLOQUES FUNCIONALES
# =============================================================

bloques <- tribble(
  ~gene,      ~bloque,
  "INTS2",    "Complejo Integrador",
  "INTS8",    "Complejo Integrador",
  "MED12",    "Complejo Mediador",
  "MED19",    "Complejo Mediador",
  "MED30",    "Complejo Mediador",
  "SUPT5H",   "Elongación transcripcional",
  "SUPT6H",   "Elongación transcripcional",
  "PAF1",     "Elongación transcripcional",
  "CTR9",     "Elongación transcripcional",
  "POLR2B",   "Elongación transcripcional",
  "HSPA9",    "Función mitocondrial",
  "PHB",      "Función mitocondrial",
  "PHB2",     "Función mitocondrial",
  "DNAJC19",  "Función mitocondrial",
  "TIMM23B",  "Función mitocondrial",
  "NCBP2",    "Procesamiento de RNA",
  "SRRT",     "Procesamiento de RNA",
  "DDX41",    "Procesamiento de RNA",
  "GATA1",    "Regulación/señalización",
  "GAB2",     "Regulación/señalización"
)

colores_bloques <- c(
  "Complejo Integrador"        = "#8e44ad",
  "Complejo Mediador"          = "#2980b9",
  "Elongación transcripcional" = "#27ae60",
  "Función mitocondrial"       = "#e67e22",
  "Procesamiento de RNA"       = "#c0392b",
  "Regulación/señalización"    = "#16a085"
)

cat("✅ Bloques funcionales definidos:\n")
print(count(bloques, bloque))
cat("\n")

# =============================================================
# 2. CARGAR DATOS
# =============================================================

chopchop     <- read_csv("data/chopchop/CHOPCHOP_top3_dataset.csv",      show_col_types = FALSE)
crisprverse  <- read_csv("data/crisprverse/crisprVerse_top3_dataset.csv", show_col_types = FALSE)
hamming      <- read_csv("results/tabla_hamming_top3.csv",               show_col_types = FALSE)
cobertura    <- read_csv("results/tabla_cobertura_umbral.csv",           show_col_types = FALSE)
posicion     <- read_csv("results/tabla_posicion_ranking.csv",           show_col_types = FALSE)
experimental <- read_csv("data/experimental/guias_experimentales_referencia.csv",
                         show_col_types = FALSE)

# Preparar CHOPCHOP limpio
cc <- chopchop %>%
  rename(gene = Gene, score = Efficiency, gc = `GC content (%)`) %>%
  left_join(bloques, by = "gene") %>%
  mutate(herramienta = "CHOPCHOP")

# Preparar crisprVerse limpio
cv <- crisprverse %>%
  rename(score = score_ruleset1) %>%
  mutate(gc = (nchar(gsub("[^GCgc]", "", protospacer)) / nchar(protospacer)) * 100) %>%
  left_join(bloques, by = "gene") %>%
  mutate(herramienta = "crisprVerse")

# =============================================================
# 3. FIGURAS DE SCORES — SEPARADAS (corrección error tutor)
# =============================================================

cat("🎨 Generando figuras de scores...\n")

# Orden de genes por bloque para que queden agrupados visualmente
orden_genes <- bloques %>% arrange(bloque, gene) %>% pull(gene)

# --- FIGURA A: Scores CHOPCHOP (escala Doench 2016, ~0-100) ---
p_cc_scores <- cc %>%
  mutate(gene = factor(gene, levels = orden_genes)) %>%
  ggplot(aes(x = gene, y = score, fill = bloque)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  geom_text(aes(label = round(score, 1)),
            position = position_dodge(0.8), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = colores_bloques) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "bottom",
        legend.title    = element_blank(),
        strip.background = element_rect(fill = "#f0f0f0")) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(
    title    = "CHOPCHOP — Scores de eficiencia por gen (top-3 guías)",
    subtitle = "Modelo: Doench 2016 Rule Set 2 | Escala: ~0 a 100\nNOTA: Esta escala NO es comparable con los scores de crisprVerse",
    x        = "Gen (agrupado por bloque funcional)",
    y        = "Score de eficiencia (Doench 2016)"
  )

ggsave("figures/figura_scores_chopchop.png",
       plot = p_cc_scores, width = 14, height = 7, dpi = 300)
cat("✅ Guardado: figures/figura_scores_chopchop.png\n")

# --- FIGURA B: Scores crisprVerse (Rule Set 1, escala 0-1) ---
p_cv_scores <- cv %>%
  mutate(gene = factor(gene, levels = orden_genes)) %>%
  ggplot(aes(x = gene, y = score, fill = bloque)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  geom_text(aes(label = round(score, 3)),
            position = position_dodge(0.8), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = colores_bloques) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "bottom",
        legend.title    = element_blank()) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(
    title    = "crisprVerse — Scores de eficiencia por gen (top-3 guías)",
    subtitle = "Modelo: Rule Set 1 (Azimuth) | Escala: 0 a 1\nNOTA: Esta escala NO es comparable con los scores de CHOPCHOP",
    x        = "Gen (agrupado por bloque funcional)",
    y        = "Score de eficiencia (Rule Set 1)"
  )

ggsave("figures/figura_scores_crisprverse.png",
       plot = p_cv_scores, width = 14, height = 7, dpi = 300)
cat("✅ Guardado: figures/figura_scores_crisprverse.png\n")

# =============================================================
# 4. FIGURA GC CONTENT (sí es comparable entre herramientas)
# =============================================================

gc_combinado <- bind_rows(
  cc %>% select(gene, bloque, gc, herramienta),
  cv %>% select(gene, bloque, gc, herramienta)
) %>% mutate(gene = factor(gene, levels = orden_genes))

p_gc <- ggplot(gc_combinado, aes(x = gene, y = gc, fill = herramienta)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6,
               alpha = 0.8, outlier.size = 1.5) +
  geom_hline(yintercept = c(40, 60), linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  annotate("text", x = 0.5, y = 41, label = "40%", color = "grey40", size = 3) +
  annotate("text", x = 0.5, y = 61, label = "60%", color = "grey40", size = 3) +
  scale_fill_manual(values = c("CHOPCHOP" = "#3498db", "crisprVerse" = "#e74c3c")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "bottom",
        legend.title    = element_blank()) +
  labs(
    title    = "Contenido GC de las guías top-3 por gen y herramienta",
    subtitle = "Rango recomendado: 40–60% (líneas discontinuas) | GC content sí es comparable entre herramientas",
    x        = "Gen",
    y        = "Contenido GC (%)"
  )

ggsave("figures/figura_gc_content.png",
       plot = p_gc, width = 14, height = 7, dpi = 300)
cat("✅ Guardado: figures/figura_gc_content.png\n")

# Tabla resumen GC
tabla_gc <- gc_combinado %>%
  group_by(herramienta, bloque) %>%
  summarise(gc_media   = round(mean(gc), 1),
            gc_mediana = round(median(gc), 1),
            gc_sd      = round(sd(gc), 1),
            n          = n(), .groups = "drop")

write_csv(tabla_gc, "results/tabla_gc_content.csv")
cat("✅ Guardado: results/tabla_gc_content.csv\n\n")

# =============================================================
# 5. ANÁLISIS POR BLOQUES — Hamming
# =============================================================

cat("🔄 Analizando concordancia por bloques funcionales...\n")

hamming_bloques <- hamming %>%
  left_join(bloques, by = "gene")

# Resumen estadístico por bloque y herramienta
tabla_bloques <- hamming_bloques %>%
  group_by(bloque, herramienta) %>%
  summarise(
    n_guias        = n(),
    hamming_media  = round(mean(hamming_minimo), 2),
    hamming_sd     = round(sd(hamming_minimo), 2),
    hamming_mediana = median(hamming_minimo),
    hamming_min    = min(hamming_minimo),
    hamming_max    = max(hamming_minimo),
    hamming_le3    = sum(hamming_minimo <= 3),
    pct_le3        = round(mean(hamming_minimo <= 3) * 100, 1),
    seed_identica  = sum(seed_identica, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n📊 TABLA BLOQUES FUNCIONALES:\n")
print(as.data.frame(tabla_bloques))
write_csv(tabla_bloques, "results/tabla_bloques_funcionales.csv")
cat("\n✅ Guardado: results/tabla_bloques_funcionales.csv\n\n")

# --- FIGURA: Hamming por bloque funcional ---
p_bloques_hamming <- hamming_bloques %>%
  ggplot(aes(x = bloque, y = hamming_minimo, fill = herramienta)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6,
               alpha = 0.85, outlier.size = 1.5) +
  geom_jitter(aes(color = herramienta),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 2.5, alpha = 0.7) +
  scale_fill_manual(values  = c("CHOPCHOP" = "#3498db", "crisprVerse" = "#e74c3c")) +
  scale_color_manual(values = c("CHOPCHOP" = "#2471a3", "crisprVerse" = "#c0392b")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x   = element_text(angle = 25, hjust = 1, face = "bold"),
        legend.position = "bottom",
        legend.title    = element_blank()) +
  guides(color = "none") +
  labs(
    title    = "Distancia Hamming mínima (top-3) por bloque funcional",
    subtitle = "Cada punto = una guía experimental | Boxplot = distribución por bloque",
    x        = "Bloque funcional",
    y        = "Hamming mínimo (nº mismatches con la mejor guía top-3)",
    fill     = "Herramienta"
  )

ggsave("figures/figura_bloques_hamming.png",
       plot = p_bloques_hamming, width = 13, height = 7, dpi = 300)
cat("✅ Guardado: figures/figura_bloques_hamming.png\n")

# =============================================================
# 6. ANÁLISIS POR BLOQUES — Cobertura en ranking completo
# =============================================================

cobertura_bloques <- posicion %>%
  left_join(bloques, by = "gene") %>%
  group_by(bloque, herramienta) %>%
  summarise(
    n_guias          = n(),
    exactas          = sum(!is.na(posicion_exacta)),
    pct_exactas      = round(mean(!is.na(posicion_exacta)) * 100, 1),
    match_h3_todos   = sum(!is.na(primera_pos_h3)),
    pct_h3_todos     = round(mean(!is.na(primera_pos_h3)) * 100, 1),
    .groups = "drop"
  )

cat("\n📊 COBERTURA POR BLOQUE (ranking completo):\n")
print(as.data.frame(cobertura_bloques))

# --- FIGURA: Cobertura exacta por bloque ---
p_bloques_cob <- cobertura_bloques %>%
  ggplot(aes(x = bloque, y = pct_exactas, fill = herramienta)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  geom_text(aes(label = paste0(pct_exactas, "%")),
            position = position_dodge(0.8), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("CHOPCHOP" = "#3498db", "crisprVerse" = "#e74c3c")) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x   = element_text(angle = 25, hjust = 1, face = "bold"),
        legend.position = "bottom",
        legend.title    = element_blank()) +
  labs(
    title    = "Coincidencia exacta en ranking completo por bloque funcional",
    subtitle = "% de guías experimentales encontradas exactamente en los rankings completos",
    x        = "Bloque funcional",
    y        = "% guías con coincidencia exacta",
    fill     = "Herramienta"
  )

ggsave("figures/figura_bloques_cobertura.png",
       plot = p_bloques_cob, width = 13, height = 7, dpi = 300)
cat("✅ Guardado: figures/figura_bloques_cobertura.png\n\n")

# =============================================================
# 7. RESUMEN FINAL EN CONSOLA
# =============================================================

cat("============================================================\n")
cat("RESUMEN FINAL — ANÁLISIS POR BLOQUES FUNCIONALES\n")
cat("============================================================\n\n")

cat("GC CONTENT — ¿Diferencias entre herramientas?\n")
tabla_gc_tool <- gc_combinado %>%
  group_by(herramienta) %>%
  summarise(media = round(mean(gc), 1), sd = round(sd(gc), 1), .groups = "drop")
print(as.data.frame(tabla_gc_tool))
cat("\n")

cat("HAMMING — ¿Qué bloque tiene mayor similitud con experimentales?\n")
tabla_bloques %>%
  group_by(bloque) %>%
  summarise(hamming_medio_global = round(mean(hamming_media), 2), .groups = "drop") %>%
  arrange(hamming_medio_global) %>%
  print()

cat("\nCOBERTURA EXACTA — ¿Qué bloque tiene más guías en los rankings?\n")
cobertura_bloques %>%
  select(bloque, herramienta, pct_exactas) %>%
  pivot_wider(names_from = herramienta, values_from = pct_exactas) %>%
  print()

cat("\n✅ Script 07 completado. Todos los archivos generados.\n")
