# =========================
# LIBRERÍAS
# =========================
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# =========================
# CARGAR DATASETS
# =========================

chopchop <- read_csv("CHOPCHOP_top3_dataset.csv")
crisprverse <- read_csv("crisprVerse_top3_dataset.csv")
experimental <- read_csv("guias_experimentales_referencia.csv")

# =========================
# HOMOGENEIZAR NOMBRES
# =========================

chopchop_clean <- chopchop %>%
  rename(
    gene = Gene,
    sequence = sgRNA_20nt,
    efficiency = Efficiency,
    gc_content = `GC content (%)`
  ) %>%
  mutate(tool = "CHOPCHOP")

crisprverse_clean <- crisprverse %>%
  rename(
    sequence = protospacer,
    efficiency = score_ruleset1
  ) %>%
  mutate(tool = "crisprVerse")

experimental_long <- experimental %>%
  select(gene, sgRNA_A_seq, sgRNA_B_seq) %>%
  pivot_longer(
    cols = c(sgRNA_A_seq, sgRNA_B_seq),
    names_to = "experimental_guide",
    values_to = "sequence"
  )

# =========================
# TABLA 1: RESUMEN GENERAL
# =========================

resumen_general <- bind_rows(
  chopchop_clean %>%
    summarise(
      herramienta = "CHOPCHOP",
      n_genes = n_distinct(gene),
      n_guias = n(),
      media_score = mean(efficiency, na.rm = TRUE),
      mediana_score = median(efficiency, na.rm = TRUE),
      sd_score = sd(efficiency, na.rm = TRUE),
      min_score = min(efficiency, na.rm = TRUE),
      max_score = max(efficiency, na.rm = TRUE)
    ),
  crisprverse_clean %>%
    summarise(
      herramienta = "crisprVerse",
      n_genes = n_distinct(gene),
      n_guias = n(),
      media_score = mean(efficiency, na.rm = TRUE),
      mediana_score = median(efficiency, na.rm = TRUE),
      sd_score = sd(efficiency, na.rm = TRUE),
      min_score = min(efficiency, na.rm = TRUE),
      max_score = max(efficiency, na.rm = TRUE)
    )
)

write.csv(resumen_general, "tabla_resumen_scores.csv", row.names = FALSE)

print(resumen_general)

# =========================
# TABLA 2: GUÍAS POR GEN
# =========================

guias_por_gen <- bind_rows(
  chopchop_clean %>%
    count(tool, gene, name = "n_guias"),
  crisprverse_clean %>%
    count(tool, gene, name = "n_guias")
)

write.csv(guias_por_gen, "tabla_guias_por_gen.csv", row.names = FALSE)

print(guias_por_gen)

# =========================
# TABLA 3: CONCORDANCIA EXACTA
# =========================

conc_chopchop <- chopchop_clean %>%
  inner_join(experimental_long, by = c("gene", "sequence"))

conc_crisprverse <- crisprverse_clean %>%
  inner_join(experimental_long, by = c("gene", "sequence"))

tabla_concordancia <- data.frame(
  herramienta = c("CHOPCHOP", "crisprVerse"),
  n_guias_top3 = c(nrow(chopchop_clean), nrow(crisprverse_clean)),
  n_guias_experimentales = c(nrow(experimental_long), nrow(experimental_long)),
  coincidencias_exactas = c(nrow(conc_chopchop), nrow(conc_crisprverse)),
  porcentaje_concordancia = c(
    nrow(conc_chopchop) / nrow(chopchop_clean) * 100,
    nrow(conc_crisprverse) / nrow(crisprverse_clean) * 100
  )
)

write.csv(tabla_concordancia, "tabla_concordancia_exacta.csv", row.names = FALSE)

print(tabla_concordancia)

# =========================
# FIGURA 1: SCORES POR HERRAMIENTA
# =========================

scores_plot <- bind_rows(
  chopchop_clean %>%
    select(gene, sequence, efficiency, tool),
  crisprverse_clean %>%
    select(gene, sequence, efficiency, tool)
)

ggplot(scores_plot, aes(x = tool, y = efficiency)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +
  theme_bw() +
  labs(
    title = "Distribución de scores de eficiencia por herramienta",
    x = "Herramienta",
    y = "Score de eficiencia"
  )

ggsave("figura_scores_por_herramienta.png", width = 7, height = 5, dpi = 300)

# =========================
# FIGURA 2: SCORES POR GEN
# =========================

ggplot(scores_plot, aes(x = gene, y = efficiency, fill = tool)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Distribución de scores de eficiencia por gen",
    x = "Gen",
    y = "Score de eficiencia"
  )

ggsave("figura_scores_por_gen.png", width = 10, height = 5, dpi = 300)

# =========================
# FIGURA 3: CONCORDANCIA EXACTA
# =========================

ggplot(tabla_concordancia, aes(x = herramienta, y = coincidencias_exactas)) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Coincidencias exactas con guías experimentales",
    x = "Herramienta",
    y = "Número de coincidencias exactas"
  )

ggsave("figura_concordancia_exacta.png", width = 7, height = 5, dpi = 300)
