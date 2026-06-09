# Benchmarking de CHOPCHOP y crisprVerse para diseño de guías CRISPR/Cas9

**Trabajo Fin de Máster — Máster Universitario en Bioinformática (UNIR)**  
**Autor:** Jaime Hernán Gil  
**Tutor:** Dr. Andrés González Jiménez  
**Título:** *Evaluación comparativa de CHOPCHOP y crisprVerse en diseño de guías CRISPR/Cas9*

---

## Descripción

Este repositorio contiene todos los scripts y datasets generados para el benchmarking comparativo de dos herramientas bioinformáticas de diseño de guías CRISPR/Cas9: CHOPCHOP (interfaz web) y crisprVerse (pipeline programático en R/Bioconductor).

El estudio utiliza como referencia experimental un subconjunto de 20 genes del dataset Perturb-seq de Replogle et al. (2022), evaluando la concordancia entre las guías diseñadas in silico y las guías empleadas experimentalmente.

---

## Estructura del repositorio

```
├── data/
│   ├── experimental/     # Guías experimentales extraídas de Replogle et al. 2022
│   ├── chopchop/         # Rankings completos y dataset top-3 de CHOPCHOP
│   └── crisprverse/      # Rankings completos y dataset top-3 de crisprVerse
├── scripts/
│   ├── 01_preprocessing.R          # Extracción y procesamiento del dataset experimental
│   ├── 02_chopchop_design.py       # Automatización del diseño de guías con CHOPCHOP v3
│   ├── 03_crisprverse_design.R     # Pipeline de diseño con crisprVerse/crisprDesign
│   ├── 04_basic_comparison.R       # Análisis comparativo básico (concordancia exacta)
│   ├── 05_concordance_analysis.R   # Análisis de similitud: distancia Hamming y región seed
│   ├── 06_suboptimal_analysis.R    # Análisis subóptimo: posición en ranking completo
│   └── 07_functional_blocks.R      # Análisis por bloques funcionales
├── figures/              # Figuras generadas (300 DPI)
├── results/              # Tablas de resultados
└── README.md
```

---

## Genes analizados (n = 20)

| Bloque funcional | Genes |
|---|---|
| Complejo Integrador | INTS2, INTS8 |
| Complejo Mediador | MED12, MED19, MED30 |
| Elongación transcripcional | SUPT5H, SUPT6H, PAF1, CTR9, POLR2B |
| Función mitocondrial | HSPA9, PHB, PHB2, DNAJC19, TIMM23B |
| Procesamiento de RNA | NCBP2, SRRT, DDX41 |
| Regulación/señalización | GATA1, GAB2 |

---

## Dataset de referencia

Los datos experimentales provienen de:

> Replogle, J.M. et al. (2022). Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell*, 185(14), 2559-2575.

La Tabla S1 del material suplementario fue utilizada como fuente de guías CRISPRi experimentales.

---

## Herramientas utilizadas

- **CHOPCHOP v3** — https://chopchop.cbu.uib.no/
- **crisprVerse** — Hoberecht et al. (2022), Bioconductor
- **R** >= 4.1.0 con paquetes: `crisprDesign`, `dplyr`, `ggplot2`, `stringdist`, `tidyr`
- **Python** >= 3.8 con: `requests`, `pandas`

---

## Referencia

Si utilizas este código, por favor cita el TFM correspondiente.
