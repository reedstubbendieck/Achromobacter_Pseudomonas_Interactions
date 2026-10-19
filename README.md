# Siderophore-Mediated Competition Between Cystic Fibrosis Pathogens Increases Virulence Factor Production

Niladri Bhowmik<sup>1</sup>, Laura Rodriguez-Velandia<sup>2</sup>, Mercedes Pérez Pérez<sup>1</sup>, Sydney M. Morabbi<sup>1</sup>, Laura Sanchez<sup>2</sup>, and Reed M. Stubbendieck<sup>1,#</sup>

<sup>1</sup>Department of Microbiology and Molecular Genetics, Oklahoma State University, Stillwater, OK 74078
<sup>2</sup>Department of Chemistry and Biochemistry, University of California, Santa Cruz, Santa Cruz, CA 95064

## Introduction

Code to replicate the computational results and figures of our study of competition between _Achromobacter_ and _Pseudomonas aeruginosa_.

## Data

| Directory | Contents |
|---|---|
| `./rawData/16S_rRNA_gene_amplicon_studies/` | SRA run accessions and clinical metadata (FEV<sub>1</sub>, age, BioProject) |
| `./rawData/clinical_pseudomonas/` | clinical _P. aeruginosa_ isolate and patient metadata |
| `./rawData/proteomics/` | MaxQuant `proteinGroups.xlsx` |
| `./rawData/rnaseq/samples.tsv` | RNA-seq sample sheet |
| `./genomes/` | _Achromobacter_ genome assemblies, one `.fna` per strain |
| `./derivedData/` | everything the scripts write |
| `./figures/`, `./tables/` | the manuscript figures and supplemental tables |

Note: The raw RNA-seq reads are in the NCBI Sequence Read Archive and are not in this repository. Set `RNASEQ_RAW` to the directory holding the FASTQ files named in `samples.tsv` and `RNASEQ_WORK` to a scratch directory for the intermediates:

```bash
export RNASEQ_RAW=/path/to/fastq
export RNASEQ_WORK=/path/to/scratch
```

## Software

| tool | version |
|---|---|
| [antiSMASH](https://docs.antismash.secondarymetabolites.org/install/) | 8.0.4 |
| [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE) | 2 |
| [clinker](https://github.com/gamcil/clinker) | |
| [fastp](https://github.com/OpenGene/fastp) | 0.23.2 |
| [kallisto](https://pachterlab.github.io/kallisto/) | 0.52.0 |
| [MEME Suite](https://meme-suite.org/meme/) | 5.5.9 |
| [Mothur](https://mothur.org/) | |
| [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | |
| [SRA Toolkit](https://github.com/ncbi/sra-tools) | |

Python: biopython, matplotlib, numpy, openpyxl, pandas, scipy.

R 4.4.3: apeglm 1.28.0, ape, circlize, ComplexHeatmap, cowplot, DESeq2 1.46.0, egg, ggrepel, ggseqlogo, ggtext, gridExtra, matrixStats, multcompView, patchwork, scales, scatterpie, tidyverse, tximport 1.34.0.

Conda environments for antiSMASH, BiG-SCAPE and clinker are in `./yml/`.

## Replicating the analysis

Run from the project root, in this order, then knit the Rmd.

### 1. Comparative genomics

antiSMASH on every genome, BiG-SCAPE to group BGCs into GCFs, run clinker on the alc BGC. `fix_gbk_organism.py` runs between the first two so BiG-SCAPE groups regions by genome. `build_presence_absence.py` collapses the clustering to one row per strain and identifies the alcaligin and DSM26089 metallophore BGCs.

```bash
python3 ./scripts/genomics/run_antismash.py --conda-env antismash
python3 ./scripts/genomics/fix_gbk_organism.py
python3 ./scripts/genomics/run_bigscape.py --conda-env bigscape --pfam /path/to/Pfam-A.hmm
python3 ./scripts/genomics/run_clinker.py --conda-env clinker
python3 ./scripts/genomics/build_presence_absence.py
```

### 2. Alcaligin pathway scan

Percent identity of each _Achromobacter_ genome's best hit to the _B. pertussis_ Tohama I alcaligin and ornithine pathway proteins.

```bash
python3 ./scripts/genomics/scan_alcaligin_pathway_genes.py
```

### 3. RNA-seq

Reference preparation, trimming, kallisto pseudoalignment, QC summary, and promoter Fur box scans. Run with `RNASEQ_RAW` and `RNASEQ_WORK` set.

```bash
bash ./scripts/rnaseq/rnaseq_run_all.sh
```

### 4. Proteomics

`proteomics_coculture_enrichment.py` reads the proteomics data, `proteomics_volcano_data.py` builds the table the Rmd plots from, and `proteomics_iron_virulence.py` writes the iron and virulence-related protein subset.

```bash
python3 ./scripts/proteomics/proteomics_coculture_enrichment.py
python3 ./scripts/proteomics/proteomics_volcano_data.py
python3 ./scripts/proteomics/proteomics_iron_virulence.py
```

### 5. Figures and tables

```bash
Rscript -e 'rmarkdown::render("Achromobacter_Pseudomonas_Interactions.Rmd")'
```
