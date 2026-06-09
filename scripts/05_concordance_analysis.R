# =============================================================
# SCRIPT 05: ANÁLISIS DE CONCORDANCIA AMPLIADA
# Distancia Hamming y análisis de región seed entre guías
# diseñadas in silico (top-3) y guías experimentales
# =============================================================
# Outputs generados:
#   results/tabla_hamming_top3.csv        — Hamming por guía y gen
#   results/tabla_resumen_hamming.csv     — Resumen por herramienta
#   results/tabla_resumen_seed.csv        — Análisis región seed
#   figures/figura_hamming_distribucion.png
#   figures/figura_heatmap_hamming.png
# =============================================================

library(dplyr)
library(readr)
library(stringdist)
library(ggplot2)
library(tidyr)

# --- Directorio de trabajo ---
BASE <- "/Users/jaimehache/Desktop/GENÉTICA/(iii) MÁSTER BIOINFORMÁTICA/(ii) SEGUNDO SEMESTRE/TFM/GITHUB"
setwd(BASE)

cat("📂 Directorio de trabajo:", BASE, "\n\n")

# =============================================================
# 1. CARGAR DATOS
# =============================================================

chopchop     <- read_csv("data/chopchop/CHOPCHOP_top3_dataset.csv",     show_col_types = FALSE)
crisprverse  <- read_csv("data/crisprverse/crisprVerse_top3_dataset.csv", show_col_types = FALSE)
experimental <- read_csv("data/experimental/guias_experimentales_referencia.csv", show_col_types = FALSE)

cat("✅ Datos cargados:\n")
cat("   CHOPCHOP top-3:   ", nrow(chopchop), "guías\n")
cat("   crisprVerse top-3:", nrow(crisprverse), "guías\n")
cat("   Experimentales:   ", nrow(experimental), "genes\n\n")

# =============================================================
# 2. PREPARAR FORMATO
# =============================================================

# Experimental: una fila por guía (A y B)
exp_long <- experimental %>%
  select(gene, sgRNA_A_seq, sgRNA_B_seq) %>%
  pivot_longer(
    cols      = c(sgRNA_A_seq, sgRNA_B_seq),
    names_to  = "guia_id",
    values_to = "seq_experimental"
  ) %>%
  mutate(guia_id = ifelse(guia_id == "sgRNA_A_seq", "A", "B"))

# CHOPCHOP: extraer columnas relevantes
cc_clean <- chopchop %>%
  rename(gene = Gene, seq_insilico = sgRNA_20nt) %>%
  select(gene, Rank, seq_insilico)

# crisprVerse: añadir rank por gen
cv_clean <- crisprverse %>%
  rename(seq_insilico = protospacer) %>%
  select(gene, seq_insilico) %>%
  group_by(gene) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

# =============================================================
# 3. FUNCIÓN PRINCIPAL: calcular Hamming mínimo
#    Para cada guía experimental, busca cuál de las top-3
#    guías in silico es la más parecida y cuántos mismatches tiene
# =============================================================

calcular_mejor_match <- function(data_tool, exp_data, nombre_tool) {

  resultados <- list()
  genes <- unique(exp_data$gene)

  for (g in genes) {
    guias_exp  <- exp_data  %>% filter(gene == g)
    guias_tool <- data_tool %>% filter(gene == g)

    if (nrow(guias_tool) == 0) {
      cat("   ⚠️  Sin datos para gen:", g, "en", nombre_tool, "\n")
      next
    }

    for (i in seq_len(nrow(guias_exp))) {
      seq_exp <- guias_exp$seq_experimental[i]
      guia_id <- guias_exp$guia_id[i]

      # Hamming entre la guía experimental y cada una de las top-3
      hammings  <- stringdist(seq_exp, guias_tool$seq_insilico, method = "hamming")
      mejor_idx <- which.min(hammings)

      resultados[[length(resultados) + 1]] <- data.frame(
        gene               = g,
        guia_exp_id        = guia_id,
        herramienta        = nombre_tool,
        seq_experimental   = seq_exp,
        mejor_seq_insilico = guias_tool$seq_insilico[mejor_idx],
        rank_mejor_match   = guias_tool$Rank[mejor_idx],
        hamming_minimo     = hammings[mejor_idx],
        hamming_rank1      = hammings[guias_tool$Rank == 1],
        hamming_rank2      = hammings[guias_tool$Rank == 2],
        hamming_rank3      = hammings[guias_tool$Rank == 3]
      )
    }
  }

  bind_rows(resultados)
}

cat("🔄 Calculando distancias Hamming...\n")
resultados_cc <- calcular_mejor_match(cc_clean, exp_long, "CHOPCHOP")
resultados_cv <- calcular_mejor_match(cv_clean, exp_long, "crisprVerse")

resultados_hamming <- bind_rows(resultados_cc, resultados_cv)

write_csv(resultados_hamming, "results/tabla_hamming_top3.csv")
cat("✅ Guardado: results/tabla_hamming_top3.csv\n\n")

# =============================================================
# 4. TABLA RESUMEN: ¿Cuántas guías tienen Hamming ≤ 1, 2, 3?
# =============================================================

resumen_hamming <- resultados_hamming %>%
  group_by(herramienta) %>%
  summarise(
    total_guias_exp = n(),
    hamming_0       = sum(hamming_minimo == 0),
    hamming_le1     = sum(hamming_minimo <= 1),
    hamming_le2     = sum(hamming_minimo <= 2),
    hamming_le3     = sum(hamming_minimo <= 3),
    pct_le1         = round(mean(hamming_minimo <= 1) * 100, 1),
    pct_le2         = round(mean(hamming_minimo <= 2) * 100, 1),
    pct_le3         = round(mean(hamming_minimo <= 3) * 100, 1),
    hamming_media   = round(mean(hamming_minimo), 2),
    hamming_mediana = median(hamming_minimo),
    hamming_min     = min(hamming_minimo),
    hamming_max     = max(hamming_minimo)
  )

cat("📊 RESUMEN HAMMING (top-3):\n")
print(as.data.frame(resumen_hamming))
write_csv(resumen_hamming, "results/tabla_resumen_hamming.csv")
cat("✅ Guardado: results/tabla_resumen_hamming.csv\n\n")

# =============================================================
# 5. ANÁLISIS REGIÓN SEED (posiciones 9-20, las 12 nt más
#    próximas al PAM — las más críticas para especificidad)
# =============================================================

resultados_hamming <- resultados_hamming %>%
  mutate(
    seed_exp      = substr(seq_experimental,   9, 20),
    seed_insilico = substr(mejor_seq_insilico, 9, 20),
    hamming_seed  = stringdist(seed_exp, seed_insilico, method = "hamming"),
    seed_identica = (hamming_seed == 0)
  )

resumen_seed <- resultados_hamming %>%
  group_by(herramienta) %>%
  summarise(
    total          = n(),
    seed_identica  = sum(seed_identica),
    pct_seed_igual = round(mean(seed_identica) * 100, 1),
    seed_hamming_media   = round(mean(hamming_seed), 2),
    seed_hamming_mediana = median(hamming_seed)
  )

cat("📊 RESUMEN REGIÓN SEED:\n")
print(as.data.frame(resumen_seed))
write_csv(resumen_seed, "results/tabla_resumen_seed.csv")

# Actualizar tabla completa con datos seed
write_csv(resultados_hamming, "results/tabla_hamming_top3.csv")
cat("✅ Guardado: results/tabla_resumen_seed.csv\n\n")

# =============================================================
# 6. FIGURAS
# =============================================================

cat("🎨 Generando figuras...\n")

# --- FIGURA A: Histograma de distribución de Hamming ---
p_hist <- ggplot(resultados_hamming, aes(x = hamming_minimo, fill = herramienta)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.85) +
  scale_x_continuous(breaks = 0:20) +
  scale_fill_manual(values = c("CHOPCHOP" = "#3498db", "crisprVerse" = "#e74c3c")) +
  facet_wrap(~herramienta, ncol = 1) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "#f0f0f0"),
        strip.text = element_text(face = "bold")) +
  labs(
    title    = "Distribución de distancias Hamming mínimas",
    subtitle = "Entre las guías top-3 in silico y las guías experimentales (Replogle et al. 2022)",
    x        = "Distancia Hamming mínima (nº de mismatches)",
    y        = "Número de guías experimentales"
  )

ggsave("figures/figura_hamming_distribucion.png",
       plot = p_hist, width = 9, height = 7, dpi = 300)
cat("✅ Guardado: figures/figura_hamming_distribucion.png\n")

# --- FIGURA B: Heatmap gen × Hamming ---
p_heat <- ggplot(resultados_hamming,
                 aes(x = gene, y = guia_exp_id, fill = hamming_minimo)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = hamming_minimo), size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "#2ecc71", high = "#e74c3c",
                      name = "Hamming\nmínimo",
                      breaks = c(0, 5, 10, 15, 20)) +
  facet_wrap(~herramienta, ncol = 1) +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
        strip.background = element_rect(fill = "#f0f0f0"),
        strip.text   = element_text(face = "bold")) +
  labs(
    title    = "Distancia Hamming mínima por gen y guía experimental",
    subtitle = "Verde = más similar, Rojo = más diferente",
    x        = "Gen",
    y        = "Guía experimental"
  )

ggsave("figures/figura_heatmap_hamming.png",
       plot = p_heat, width = 14, height = 7, dpi = 300)
cat("✅ Guardado: figures/figura_heatmap_hamming.png\n\n")

# =============================================================
# 7. RESUMEN FINAL EN CONSOLA
# =============================================================

cat("============================================================\n")
cat("RESUMEN FINAL — ANÁLISIS DE CONCORDANCIA AMPLIADA\n")
cat("============================================================\n\n")

for (tool in c("CHOPCHOP", "crisprVerse")) {
  r <- resumen_hamming %>% filter(herramienta == tool)
  s <- resumen_seed    %>% filter(herramienta == tool)
  cat(tool, ":\n")
  cat("  Guías evaluadas:          ", r$total_guias_exp, "\n")
  cat("  Coincidencia exacta:       ", r$hamming_0, "(", r$hamming_0/r$total_guias_exp*100, "%)\n")
  cat("  Hamming ≤ 1 mismatch:     ", r$hamming_le1, "(", r$pct_le1, "%)\n")
  cat("  Hamming ≤ 2 mismatches:   ", r$hamming_le2, "(", r$pct_le2, "%)\n")
  cat("  Hamming ≤ 3 mismatches:   ", r$hamming_le3, "(", r$pct_le3, "%)\n")
  cat("  Hamming media / mediana:  ", r$hamming_media, "/", r$hamming_mediana, "\n")
  cat("  Región seed idéntica:     ", s$seed_identica, "(", s$pct_seed_igual, "%)\n\n")
}

cat("✅ Script 05 completado. Archivos generados en results/ y figures/\n")
