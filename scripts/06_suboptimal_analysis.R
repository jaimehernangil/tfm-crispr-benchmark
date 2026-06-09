# =============================================================
# SCRIPT 06: ANÁLISIS SUBÓPTIMO
# ¿En qué posición del ranking COMPLETO aparece la guía
# experimental? ¿A partir de qué top-N habría coincidencia?
# =============================================================
# Usa los rankings completos (no solo top-3):
#   - CHOPCHOP:    hasta ~500 guías por gen
#   - crisprVerse: hasta ~21.000 guías por gen
#
# Outputs:
#   results/tabla_posicion_ranking.csv    — posición exacta y Hamming por guía
#   results/tabla_cobertura_umbral.csv    — % cobertura según top-N
#   figures/figura_curva_cobertura.png
#   figures/figura_posicion_ranking.png
# =============================================================

library(dplyr)
library(readr)
library(stringdist)
library(ggplot2)
library(tidyr)
library(purrr)

BASE <- "/Users/jaimehache/Desktop/GENÉTICA/(iii) MÁSTER BIOINFORMÁTICA/(ii) SEGUNDO SEMESTRE/TFM/GITHUB"
setwd(BASE)

cat("📂 Directorio:", BASE, "\n\n")

# =============================================================
# 1. CARGAR DATOS EXPERIMENTALES
# =============================================================

experimental <- read_csv("data/experimental/guias_experimentales_referencia.csv",
                         show_col_types = FALSE)

exp_long <- experimental %>%
  select(gene, sgRNA_A_seq, sgRNA_B_seq) %>%
  pivot_longer(cols      = c(sgRNA_A_seq, sgRNA_B_seq),
               names_to  = "guia_id",
               values_to = "seq_experimental") %>%
  mutate(guia_id = ifelse(guia_id == "sgRNA_A_seq", "A", "B"))

genes <- sort(unique(exp_long$gene))
cat("Genes a analizar:", length(genes), "\n",
    paste(genes, collapse = ", "), "\n\n")

# =============================================================
# 2. CARGAR RANKINGS COMPLETOS
# =============================================================

cat("🔄 Cargando rankings completos de CHOPCHOP...\n")

leer_chopchop_completo <- function(gen) {
  archivo <- paste0("data/chopchop/", gen, "_tabla.txt")
  if (!file.exists(archivo)) { cat("  ⚠️  No encontrado:", archivo, "\n"); return(NULL) }
  df <- read_tsv(archivo, show_col_types = FALSE)
  df %>%
    mutate(gene        = gen,
           Rank        = row_number(),
           seq_insilico = substr(`Target sequence`, 1, 20)) %>%
    select(gene, Rank, seq_insilico)
}

chopchop_all <- map_dfr(genes, leer_chopchop_completo)
cat("  ✅ CHOPCHOP: ", nrow(chopchop_all), "guías en total\n")
cat("  Rango por gen:", chopchop_all %>% count(gene) %>% summarise(min=min(n), max=max(n)) %>% paste(collapse=" - "), "\n\n")

cat("🔄 Cargando rankings completos de crisprVerse...\n")

leer_crisprverse_completo <- function(gen) {
  archivo <- paste0("data/crisprverse/", gen, "_crisprVerse_ALL.csv")
  if (!file.exists(archivo)) { cat("  ⚠️  No encontrado:", archivo, "\n"); return(NULL) }
  df <- read_csv(archivo, show_col_types = FALSE)
  df %>%
    arrange(desc(score_ruleset1)) %>%
    mutate(gene         = gen,
           Rank         = row_number(),
           seq_insilico = protospacer) %>%
    select(gene, Rank, seq_insilico)
}

crisprverse_all <- map_dfr(genes, leer_crisprverse_completo)
cat("  ✅ crisprVerse:", nrow(crisprverse_all), "guías en total\n\n")

# =============================================================
# 3. FUNCIÓN PRINCIPAL: buscar posición en ranking completo
# =============================================================

analizar_posicion <- function(ranking_all, exp_data, nombre_tool) {

  resultados <- list()
  cat("  Procesando", nombre_tool, "...\n")

  for (g in genes) {
    guias_exp <- exp_data    %>% filter(gene == g)
    ranking_g <- ranking_all %>% filter(gene == g)

    if (nrow(ranking_g) == 0) next

    seqs_ranking <- ranking_g$seq_insilico

    for (i in seq_len(nrow(guias_exp))) {
      seq_exp <- guias_exp$seq_experimental[i]
      guia_id <- guias_exp$guia_id[i]

      # Coincidencia exacta en el ranking completo
      pos_exacta <- which(seqs_ranking == seq_exp)
      posicion_exacta <- ifelse(length(pos_exacta) > 0, min(pos_exacta), NA_integer_)

      # Distancia Hamming contra todo el ranking
      hammings <- stringdist(seq_exp, seqs_ranking, method = "hamming")

      mejor_hamming <- min(hammings)
      pos_mejor     <- which.min(hammings)

      # Primera posición con cada umbral de Hamming
      primera_h1 <- { p <- which(hammings <= 1); ifelse(length(p) > 0, min(p), NA_integer_) }
      primera_h2 <- { p <- which(hammings <= 2); ifelse(length(p) > 0, min(p), NA_integer_) }
      primera_h3 <- { p <- which(hammings <= 3); ifelse(length(p) > 0, min(p), NA_integer_) }
      primera_h5 <- { p <- which(hammings <= 5); ifelse(length(p) > 0, min(p), NA_integer_) }

      resultados[[length(resultados) + 1]] <- data.frame(
        gene              = g,
        guia_id           = guia_id,
        herramienta       = nombre_tool,
        total_candidatos  = nrow(ranking_g),
        posicion_exacta   = posicion_exacta,
        mejor_hamming     = mejor_hamming,
        posicion_mejor    = pos_mejor,
        primera_pos_h1    = primera_h1,
        primera_pos_h2    = primera_h2,
        primera_pos_h3    = primera_h3,
        primera_pos_h5    = primera_h5
      )
    }
  }

  bind_rows(resultados)
}

cat("🔄 Analizando posiciones en ranking completo...\n")
cat("   (crisprVerse puede tardar 1-2 minutos por el volumen de datos)\n\n")

res_cc <- analizar_posicion(chopchop_all,    exp_long, "CHOPCHOP")
res_cv <- analizar_posicion(crisprverse_all, exp_long, "crisprVerse")

resultados_posicion <- bind_rows(res_cc, res_cv)
write_csv(resultados_posicion, "results/tabla_posicion_ranking.csv")
cat("✅ Guardado: results/tabla_posicion_ranking.csv\n\n")

# =============================================================
# 4. RESUMEN: coincidencias en el ranking completo
# =============================================================

resumen_posicion <- resultados_posicion %>%
  group_by(herramienta) %>%
  summarise(
    total_guias_exp     = n(),
    total_candidatos_media = round(mean(total_candidatos)),
    coincidencia_exacta = sum(!is.na(posicion_exacta)),
    mejor_hamming_media = round(mean(mejor_hamming), 2),
    mejor_hamming_min   = min(mejor_hamming),
    match_h1_existe     = sum(!is.na(primera_pos_h1)),
    match_h2_existe     = sum(!is.na(primera_pos_h2)),
    match_h3_existe     = sum(!is.na(primera_pos_h3)),
    match_h5_existe     = sum(!is.na(primera_pos_h5))
  )

cat("📊 RESUMEN — RANKING COMPLETO:\n")
print(as.data.frame(resumen_posicion))
cat("\n")

# =============================================================
# 5. ANÁLISIS DE COBERTURA POR UMBRAL TOP-N
#    Para cada umbral top-N, qué % de guías experimentales
#    tienen al menos un match con Hamming ≤ 3 en las top-N
# =============================================================

umbrales <- c(3, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000)

tabla_umbrales <- map_dfr(umbrales, function(u) {
  resultados_posicion %>%
    mutate(match = !is.na(primera_pos_h3) & primera_pos_h3 <= u) %>%
    group_by(herramienta) %>%
    summarise(
      umbral_topN     = u,
      n_guias_exp     = n(),
      guias_con_match = sum(match),
      porcentaje      = round(mean(match) * 100, 1),
      .groups = "drop"
    )
})

# Añadir fila "TODOS"
fila_todos <- resultados_posicion %>%
  mutate(match = !is.na(primera_pos_h3)) %>%
  group_by(herramienta) %>%
  summarise(umbral_topN = max(total_candidatos),
            n_guias_exp = n(),
            guias_con_match = sum(match),
            porcentaje = round(mean(match) * 100, 1),
            .groups = "drop")

tabla_umbrales <- bind_rows(tabla_umbrales, fila_todos) %>%
  arrange(herramienta, umbral_topN)

cat("📊 COBERTURA POR UMBRAL (Hamming ≤ 3):\n")
print(as.data.frame(tabla_umbrales))
write_csv(tabla_umbrales, "results/tabla_cobertura_umbral.csv")
cat("\n✅ Guardado: results/tabla_cobertura_umbral.csv\n\n")

# =============================================================
# 6. FIGURAS
# =============================================================

cat("🎨 Generando figuras...\n")

# --- FIGURA A: Curva de cobertura según top-N ---
umbrales_plot <- tabla_umbrales %>%
  filter(umbral_topN <= 10000)

p_curva <- ggplot(umbrales_plot,
                  aes(x = umbral_topN, y = porcentaje, color = herramienta)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = c(3, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000),
    labels = c("top-3","top-5","top-10","top-20","top-50",
               "top-100","top-200","top-500","top-1000","top-5000","top-10000")
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_color_manual(values = c("CHOPCHOP" = "#3498db", "crisprVerse" = "#e74c3c")) +
  theme_bw(base_size = 13) +
  theme(axis.text.x  = element_text(angle = 40, hjust = 1),
        legend.position = "bottom",
        legend.title    = element_blank()) +
  labs(
    title    = "Cobertura acumulada de guías experimentales según umbral top-N",
    subtitle = "Criterio: al menos un match con Hamming ≤ 3 mismatches en las top-N candidatas",
    x        = "Umbral top-N (escala logarítmica)",
    y        = "% guías experimentales cubiertas"
  )

ggsave("figures/figura_curva_cobertura.png",
       plot = p_curva, width = 11, height = 7, dpi = 300)
cat("✅ Guardado: figures/figura_curva_cobertura.png\n")

# --- FIGURA B: Mejor Hamming en ranking completo, por gen ---
p_mejor <- resultados_posicion %>%
  ggplot(aes(x = reorder(gene, mejor_hamming), y = mejor_hamming,
             fill = herramienta)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "#27ae60", linewidth = 0.7) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "#e67e22", linewidth = 0.7) +
  annotate("text", x = 0.6, y = 3.3, label = "Hamming = 3", color = "#27ae60", size = 3.5) +
  annotate("text", x = 0.6, y = 5.3, label = "Hamming = 5", color = "#e67e22", size = 3.5) +
  scale_fill_manual(values = c("CHOPCHOP" = "#3498db", "crisprVerse" = "#e74c3c")) +
  facet_wrap(~guia_id, ncol = 1, labeller = labeller(guia_id = c(A = "Guía A", B = "Guía B"))) +
  theme_bw(base_size = 12) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.background = element_rect(fill = "#f0f0f0"),
        strip.text = element_text(face = "bold")) +
  labs(
    title    = "Mejor distancia Hamming en el ranking completo por gen",
    subtitle = "Distancia mínima encontrada en todos los candidatos generados",
    x        = "Gen (ordenado por Hamming medio)",
    y        = "Distancia Hamming mínima",
    fill     = "Herramienta"
  )

ggsave("figures/figura_mejor_hamming_ranking.png",
       plot = p_mejor, width = 14, height = 8, dpi = 300)
cat("✅ Guardado: figures/figura_mejor_hamming_ranking.png\n\n")

# =============================================================
# 7. RESUMEN FINAL EN CONSOLA
# =============================================================

cat("============================================================\n")
cat("RESUMEN FINAL — ANÁLISIS SUBÓPTIMO\n")
cat("============================================================\n\n")

for (tool in c("CHOPCHOP", "crisprVerse")) {
  r  <- resumen_posicion %>% filter(herramienta == tool)
  u3 <- tabla_umbrales   %>% filter(herramienta == tool, umbral_topN == 3)
  u10  <- tabla_umbrales %>% filter(herramienta == tool, umbral_topN == 10)
  u100 <- tabla_umbrales %>% filter(herramienta == tool, umbral_topN == 100)
  utodos <- fila_todos   %>% filter(herramienta == tool)

  cat(tool, "(media candidatos por gen:", r$total_candidatos_media, ")\n")
  cat("  Coincidencia exacta en ranking completo:", r$coincidencia_exacta, "\n")
  cat("  Mejor Hamming encontrado (mínimo):", r$mejor_hamming_min, "\n")
  cat("  Mejor Hamming medio:", r$mejor_hamming_media, "\n")
  cat("  Con Hamming ≤ 3 en top-3:  ", u3$guias_con_match, "/", u3$n_guias_exp, "(", u3$porcentaje, "%)\n")
  cat("  Con Hamming ≤ 3 en top-10: ", u10$guias_con_match, "/", u10$n_guias_exp, "(", u10$porcentaje, "%)\n")
  cat("  Con Hamming ≤ 3 en top-100:", u100$guias_con_match, "/", u100$n_guias_exp, "(", u100$porcentaje, "%)\n")
  cat("  Con Hamming ≤ 3 en TODOS:  ", utodos$guias_con_match, "/", utodos$n_guias_exp, "(", utodos$porcentaje, "%)\n\n")
}

cat("✅ Script 06 completado.\n")
