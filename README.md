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

## Estructura del proyecto

```
balderas-ncrnas-hypoxia-ipf/
├── README.md
├── LICENSE
├── .gitignore
├── balderas-ncrnas-hypoxia-ipf.Rproj
├── scripts/                    # Scripts de análisis de redes ncRNA
│   ├── encori.R                # Consultas a ENCORI/StarBase
│   ├── encori2.R
│   ├── network.R               # Construcción de redes lncRNA-mRNA
│   ├── rna-rna-function.R      # Funciones para interacciones RNA-RNA
│   └── run-rna-rna-from-list.R
├── data/                       # Tablas de interacciones RNA
├── results/                    # Redes, figuras
└── docs/                       # Manuscrito en preparación
```

## Descripción

A partir de los genes diferencialmente expresados identificados en el análisis transcriptómico corregido (Clariom D, hipoxia vs normoxia en fibroblastos IPF y control), se realizó un análisis de interacciones RNA-RNA utilizando ENCORI/StarBase para construir redes de regulación lncRNA-mRNA. Las redes se visualizaron y analizaron en Cytoscape.

## Herramientas

- R / RStudio
- ENCORI / StarBase (interacciones RNA-RNA)
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
