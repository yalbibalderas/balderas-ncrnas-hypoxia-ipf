# balderas-ncrnas-hypoxia-ipf

Análisis de redes de RNA no codificante (lncRNA-mRNA) en fibroblastos pulmonares de pacientes con fibrosis pulmonar idiopática (IPF) bajo condiciones de hipoxia y normoxia. Basado en datos de microarreglos Clariom D (Affymetrix) del estudio transcriptómico publicado en Cells 2022.

> **Datos base:** Los datos de expresión provienen del análisis corregido (v2) del repositorio [romero-balderas-hypoxia-ipf-fibroblasts](https://github.com/yalbibalderas/romero-balderas-hypoxia-ipf-fibroblasts). Los archivos CEL crudos están disponibles en FigShare (ver referencia en el repo original).

## Estado

🔬 **Paper en preparación**

## Colaboradores

- Dra. Yalbi I. Balderas-Martínez (Laboratorio de Biología Computacional, INER)
- Dr. Iván Salido (Instituto Nacional de Perinatología)
- Dr. Yair Romero (Facultad de Ciencias, UNAM)
- Dr. Arnoldo Aquino-Gálvez (Departamento de Fibrosis Pulmonar, INER)

## Versiones del análisis

Este repositorio contiene dos versiones del análisis, organizadas por la procedencia de los datos de expresión diferencial:

- **v1-original-2022/** — Análisis original realizado por el Dr. Iván Salido usando `affy::rma()` con diseño experimental hardcoded. Las consultas ENCORI y redes MirNet se construyeron sobre estas tablas DE.
- **v2-corrected-2025/** — Análisis actualizado usando las tablas DE corregidas generadas con `oligo::rma()` y metadata.csv (script 08 del repo de transcriptómica). Incluye 5 comparaciones DE y figuras actualizadas.

> Para detalles sobre las diferencias entre versiones, ver [CHANGELOG.md](CHANGELOG.md).

## Estructura del proyecto

```
balderas-ncrnas-hypoxia-ipf/
├── README.md, CHANGELOG.md, LICENSE, .gitignore
├── balderas-ncrnas-hypoxia-ipf.Rproj
│
├── v1-original-2022/                    # Análisis original (affy::rma)
│   ├── scripts/
│   │   ├── ivan/
│   │   │   ├── 00-pipeline-ipf-hipoxia-ivan.R   # Pipeline monolítico DE + ncRNA
│   │   │   └── rhistory-*.R                      # Historiales de sesiones R
│   │   ├── encori.R                     # Consultas ENCORI/StarBase (v1)
│   │   ├── network.R                    # Integración con RISE database
│   │   └── rna-rna-function.R           # Función rna_rna() con URLencode
│   ├── data/
│   │   ├── *-table-differential-expression-nr.csv  # Tablas DE Iván (IPF y Control)
│   │   ├── rna-interactions-resultados-*.csv       # Resultados ENCORI
│   │   ├── deg-lists/                   # Listas de DEGs filtrados
│   │   │   ├── fibrotic-hx-vs-nx/
│   │   │   └── normal-hx-vs-nx/
│   │   ├── mirnet/                      # Datos de redes MirNet
│   │   └── lncrna-rnacentral/          # Anotaciones RNAcentral
│   ├── results/
│   │   ├── enrichment/normal/mirna/     # Enrichr por miRNA individual
│   │   ├── figures/
│   │   │   ├── fig1-experimental-design.png
│   │   │   ├── heatmaps-volcanos/       # Heatmaps y volcanos (v1)
│   │   │   ├── networks/               # Screenshots de redes
│   │   │   ├── venns/                   # Diagramas de Venn
│   │   │   └── previous-versions/       # Versiones anteriores
│   │   ├── networks/                    # Redes Cytoscape (.cys)
│   │   └── tables/                      # Tablas suplementarias
│   └── docs/
│       ├── manuscript/                  # Borradores del manuscrito
│       └── photo-*.jpg                  # Fotos de reuniones
│
└── v2-corrected-2025/                   # Análisis corregido (oligo::rma)
    ├── scripts/
    │   ├── encori2.R                    # ENCORI mejorado (httr, error handling)
    │   └── run-rna-rna-from-list.R      # Consultas ENCORI desde listas DE
    ├── data/
    │   ├── deg-tables/                  # 5 tablas DE corregidas (NR_ ncRNA)
    │   ├── expression-matrices/         # Matrices de expresión normalizadas
    │   └── encori-reference/            # Tablas de referencia ENCORI/StarBase
    ├── results/
    │   └── figures/
    │       ├── fig3-heatmap-volcano-boxplot.png  # Figura compuesta (Iván)
    │       ├── heatmaps-volcanos/       # Heatmaps y volcanos corregidos
    │       └── networks/                # Figuras de redes actualizadas
    └── docs/
        └── legend-figure-2-ivan.pdf     # Leyenda Fig. 2 (Iván)
```

## Descripción

A partir de los genes diferencialmente expresados identificados en el análisis transcriptómico (Clariom D, hipoxia vs normoxia en fibroblastos IPF y control), se realizó un análisis de interacciones RNA-RNA utilizando ENCORI/StarBase para construir redes de regulación lncRNA-mRNA. Las redes se visualizaron y analizaron en Cytoscape.

### Flujo de análisis

1. **Expresión diferencial** — Identificación de ncRNAs (accession NR_) diferencialmente expresados entre condiciones de hipoxia y normoxia
2. **Interacciones RNA-RNA** — Consultas a ENCORI/StarBase para identificar interacciones lncRNA-mRNA y miRNA-target
3. **Redes MirNet** — Construcción de redes miRNA-target usando MirNet
4. **Enriquecimiento funcional** — Análisis de enriquecimiento por miRNA individual usando Enrichr (GO Biological Process)
5. **Visualización** — Redes de interacción en Cytoscape, heatmaps, volcanos y diagramas de Venn

## Herramientas

- R / RStudio
- ENCORI / StarBase (interacciones RNA-RNA)
- MirNet (redes miRNA-target)
- Enrichr (análisis de enriquecimiento funcional)
- Cytoscape (visualización de redes)
- Clariom D (Affymetrix) — datos base del análisis transcriptómico

## Datos relacionados

- Análisis transcriptómico (datos base): [romero-balderas-hypoxia-ipf-fibroblasts](https://github.com/yalbibalderas/romero-balderas-hypoxia-ipf-fibroblasts)
- Proteómica mitocondrial (mismo sistema): [romero-mitochondrial-proteomics-ipf](https://github.com/yalbibalderas/romero-mitochondrial-proteomics-ipf)

---

## English

Non-coding RNA network analysis (lncRNA-mRNA) in lung fibroblasts from idiopathic pulmonary fibrosis (IPF) patients under hypoxia and normoxia conditions. Based on Clariom D microarray data from the transcriptomic study published in Cells 2022.

**Status:** Paper in preparation.

### Related repositories

- Transcriptomic analysis (base data): [romero-balderas-hypoxia-ipf-fibroblasts](https://github.com/yalbibalderas/romero-balderas-hypoxia-ipf-fibroblasts)
- Mitochondrial proteomics (same system): [romero-mitochondrial-proteomics-ipf](https://github.com/yalbibalderas/romero-mitochondrial-proteomics-ipf)
