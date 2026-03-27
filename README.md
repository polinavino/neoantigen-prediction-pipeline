# Neoantigen Prediction Pipeline

Reproduction of a personalized cancer vaccine neoantigen prediction pipeline 
using TCGA melanoma (TCGA-SKCM) whole exome sequencing data.

## What this does

Takes somatic mutation data from a tumor sample and predicts which mutations 
produce peptides (neoantigens) likely to bind MHC class I molecules and trigger 
an immune response — the basis of personalized cancer vaccines.

## Pipeline steps

1. Download somatic mutation MAF from GDC (TCGA-SKCM project)
2. Convert MAF to VCF format (`maf_to_vep_vcf.py`)
3. Annotate VCF with Ensembl VEP including Wildtype and Frameshift plugins
4. Run pVACseq with NetMHCpan binding predictions via IEDB REST API

## Tools required

- Python 3.10
- pVACtools 5.3.0
- Ensembl VEP 115
- BWA-MEM, samtools (for alignment if starting from FASTQ)
- Conda (environment management)

## Data

Input data: TCGA-SKCM masked somatic mutation MAF, downloaded from 
https://portal.gdc.cancer.gov/projects/TCGA-SKCM
Patient sample: TCGA-ER-A2NC (file ID: 34745fc1-1ae6-4311-a7a3-fd43eeb40e4e)

Note: HLA type was assumed as HLA-A*02:01 for demonstration. In a real 
clinical application, HLA type would be determined from sequencing data 
using OptiType or similar.

## Results

48 candidate neoantigens predicted. Top candidates by fold change 
(mutation-induced improvement in MHC binding):

| Gene | Mutant Peptide | IC50 (nM) | Fold Change |
|------|---------------|-----------|-------------|
| ALDH8A1 | ILSDPLVSI | 12.36 | 1355 |
| FAM151A | SLGWTTFYM | 22.02 | 670 |
| CNGB3 | KLLTKVNKV | 21.15 | 302 |
| CDH7 | KMSPVGTSV | 21.46 | 75 |
| SLC6A3 | MAMVPIYAV | 27.83 | 7.8 |

## Known limitations

- HLA type assumed rather than typed from sequencing
- Only MHC class I (CD8+ T cell) epitopes predicted; class II not included
- Single HLA allele tested; real patients have 6 class I alleles
- NetMHCpan predictions via IEDB REST API; local IEDB install would be faster
- No RNA expression filtering applied (expression data not available in this dataset)

## Installation notes

Several patches were required to run pVACtools 5.3.0 on Apple Silicon macOS:
- Fixed asparagine/aspartate naming mismatch in manufacturability calculations
- Updated IEDB API URLs from HTTP to HTTPS
- Added request timeouts to prevent hanging on slow API responses
- Fixed Python interpreter path for subprocess calls
