# Changelog

## Diferencias entre versiones / Version differences

### v1-original-2022

Análisis original realizado entre 2022-2023 por el Dr. Iván Salido.

**Normalización:**
- `affy::rma()` con diseño experimental hardcoded: `rep(c("IPF", "Ctr"), each = 18)`
- Ruta de datos hardcoded: `/Volumes/GoogleDrive/My Drive/CEL_ARNOLDO_INER`

**Expresión diferencial:**
- 2 contrastes: Fibrotic Hx vs Nx, Normal (Control) Hx vs Nx
- Pipeline monolítico con `difExp.limma()` (función personalizada)
- Filtros: adj.P.Val < 0.05, |logFC| > 0.7 para listas de ncRNAs

**Interacciones RNA-RNA:**
- Consultas ENCORI/StarBase usando `utils::read.csv(url(link))`
- Interacciones lncRNA-mRNA y miRNA-target
- Resultados almacenados en CSVs planos

**Redes y enriquecimiento:**
- MirNet para redes miRNA-target (nodos y targets)
- Enrichr para GO Biological Process por miRNA individual
- Cytoscape para visualización de redes (.cys)

### v2-corrected-2025

Análisis actualizado usando las tablas DE corregidas del repositorio de transcriptómica.

**Normalización (heredada del repo de transcriptómica):**
- `oligo::rma()` con metadata desde `data/metadata.csv`
- Diseño experimental construido programáticamente desde metadata
- Corrección del error de asignación de muestras a condiciones

**Expresión diferencial:**
- 5 contrastes completos (todas las combinaciones relevantes):
  1. Fibrotic Hipoxia vs Fibrotic Normoxia
  2. Normal Hipoxia vs Normal Normoxia
  3. Fibrotic Hipoxia vs Normal Hipoxia
  4. Fibrotic Hipoxia vs Normal Normoxia
  5. Fibrotic Normoxia vs Normal Normoxia
- Filtrado por accession NR_ (NCBI non-coding RNA)

**Interacciones RNA-RNA:**
- Script mejorado `encori2.R` con `httr::GET()`, timeouts y manejo de errores
- `run-rna-rna-from-list.R` para consultas batch desde listas DE

**Figuras:**
- Heatmaps y volcanos regenerados con datos corregidos
- Figuras de redes actualizadas (Fig3 panels A-C)

### Diferencias clave

| Aspecto | v1-original | v2-corrected |
|---------|-------------|--------------|
| Normalización | `affy::rma()` | `oligo::rma()` |
| Diseño experimental | Hardcoded | Desde metadata.csv |
| Contrastes DE | 2 | 5 |
| Consultas ENCORI | `read.csv(url())` | `httr::GET()` con error handling |
| Matrices de expresión | Idénticas | Idénticas |
| Tablas DE | Difieren por error en asignación | Corregidas |

> **Nota:** Las matrices de expresión normalizadas son idénticas entre v1 y v2 (verificado por MD5). La diferencia surge en la asignación de muestras a condiciones para el análisis de expresión diferencial.
