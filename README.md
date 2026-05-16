# Hybrid_Rattlesnakes

Bioinformatics pipelines and analysis scripts for trio-based detection of Mendelian inheritance errors (MIEs) and *de novo* mutations (DNMs) in F1 *Crotalus viridis* × *C. scutulatus* hybrid offspring.

This repository accompanies a manuscript currently in preparation. Sequencing data accessions and a citation will be added once the paper is submitted.

---

## About

We performed whole-genome resequencing on a controlled experimental cross between a *C. viridis* dam and a *C. scutulatus* sire (the most deeply divergent species pair in the Western Rattlesnake species complex) and their three F1 hybrid offspring. The scripts in this repository implement the full bioinformatics pipeline used to quantify Mendelian inheritance patterns and to identify *de novo* mutations across the genome.

The pipeline is built around two distinct analytical tracks that share an upstream read-processing and alignment workflow:

- **Mendelian inheritance error detection** uses GATK's GVCF-mode joint genotyping followed by per-trio extraction and `bcftools +mendelian2`.
- ***De novo* mutation detection** uses multi-sample trio calling with GATK HaplotypeCaller, followed by the GATK Genotype Refinement workflow with `PossibleDeNovo` annotation and a read-level filtering step adapted from Zhang *et al.* (2020).

Scripts were originally developed and run on a SLURM-managed HPC cluster, but the cleaned versions in this repository are written as reusable templates with configurable paths, sample lists, and filtering thresholds.

---

## Pipeline overview

The scripts are numbered to reflect the order in which they were run:

| Step | Script | Description |
|------|--------|-------------|
| 0 | `00_prepare_reference.sh` | Integrate the W chromosome (Schield *et al.* 2022) into the CroVir_3.0 reference (Schield *et al.* 2019) and index. One-time setup. |
| 1a | `01a_trim_offspring_dam.sh` | Trimmomatic on paired-end reads from the offspring and dam. |
| 1b | `01b_trim_sire.sh` | Trimmomatic on multi-lane paired-end reads from the sire (received separately from a different sequencing facility; parameters differ from 1a). |
| 2 | `02_bwa_align.sh` | BWA-MEM alignment, sort, multi-lane merge, Picard MarkDuplicates. |
| 3a | `03a_repeatmasker.sh` | Build a repeat-region mask BED with RepeatMasker. |
| 3b | `03b_cnvnator.sh` | Per-sample CNV calling with CNVnator. Outputs per-sample BEDs; merge across samples downstream. |
| 4 | `04_mie_pipeline.sh` | Full MIE pipeline: GVCF-mode joint genotyping, hard filtering, masking, per-trio MIE detection. |
| 5a | `05a_dnm_trio_call.sh` | Multi-sample trio HaplotypeCaller for each offspring. |
| 5b | `05b_dnm_refine.sh` | Hard filtering, GATK Genotype Refinement, `PossibleDeNovo` annotation. |
| 5c | `05c_dnm_readlevel_filter.sh` | Read-level DNM filtering via `filter_dnm_candidates.py`. |
| 6a | `06a_snpeff_build.sh` | Build a custom SnpEff database from the reference FASTA and a MAKER GFF. One-time setup. |
| 6b | `06b_snpeff_annotate.sh` | Annotate final DNMs with SnpEff; uses `dnm_tsv_to_vcf.py` and `parse_snpeff_results.py`. |
| 7a | `07a_mie_circular_manhattan.R` | Circular Manhattan plot of MIE density (350 kb windows). |
| 7b | `07b_mie_callable_correlation.R` | MIE count vs. callable sites — correlation, regression, and triangle-layout scatter plot. |

### Helper scripts

**Python:**

- `filter_dnm_candidates.py` — read-level filter for high-confidence *de novo* candidates (DP, parental allele, allele balance).
- `dnm_tsv_to_vcf.py` — converts the final DNM TSV to a minimal VCF for SnpEff input.
- `parse_snpeff_results.py` — parses SnpEff-annotated VCFs into a human-readable summary TSV.

---

## Requirements

Tested with the versions listed; other recent versions are likely to work but have not been verified.

| Tool | Version |
|------|---------|
| Trimmomatic | 0.39 |
| BWA | 0.7.18 |
| samtools / bcftools | 1.20 |
| Picard | (any recent) |
| RepeatMasker | 4.1.9 |
| CNVnator | 0.4.1 |
| GATK | 4.6.1.0 |
| vcftools | (any recent) |
| bedtools | (any recent) |
| SnpEff | 5.4a |
| Python | ≥ 3.8 |
| R | ≥ 4.0 |

**R packages used by the analysis scripts:** `tidyverse`, `circlize`, `patchwork`, `scales`. The scripts install these on first run if missing.

A SLURM scheduler is assumed for batch submission of the bash scripts, but they will also run in a plain shell if you set `SLURM_CPUS_PER_TASK` manually and submit each script directly with `bash`.

---

## Usage

Each script is self-documenting in its header comment. The general pattern is:

1. Open the script and edit the configuration block at the top — paths, sample IDs, and any filtering thresholds you want to override.
2. Submit with `sbatch <script>.sh` (or, for scripts using a SLURM array, adjust the `--array=0-N` line to match the number of samples or offspring). The R scripts run interactively or with `Rscript`.
3. Wait for one step to complete before submitting the next; the scripts assume the outputs of the preceding step are in place.

The pipeline order is:

```
00 → 01a / 01b → 02 → 03a / 03b → 04 → 07a, 07b    (MIE track + figures)
                                  → 05a → 05b → 05c → 06a → 06b   (DNM track)
```

Steps 03a, 03b, and 04 are independent of the DNM track and can be run in parallel with it once alignment (02) is complete. The R analysis scripts (07a, 07b) run on the outputs of step 4 and are independent of the DNM track.

---

## Notes on the data

A few project-specific details worth knowing if you want to reproduce or adapt this pipeline:

- **Two sequencing batches.** The dam and three offspring were sequenced in one batch; the sire was sequenced separately at a different facility. The two trimming scripts (`01a` and `01b`) reflect the parameter choices made for each batch.
- **Reference genome.** The CroVir_3.0 reference (Schield *et al.* 2019) does not include the W chromosome; we integrated the W assembly from Schield *et al.* 2022 to support analyses of the female samples (one offspring is female). Script `00_prepare_reference.sh` documents how.
- **Coverage.** Sequencing depth was approximately 15× per sample.

---

## Data availability

- **Raw reads:** [SRA BioProject accession to be added upon submission]
- **Reference genome:** publicly available from NCBI ([GCA_003400415.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003400415.2/), Schield *et al.* 2019)
- **W chromosome assembly:** publicly available from NCBI ([GCA_024760675.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_024760675.1/), Schield *et al.* 2022)

---

## Citation

If you use this code, please cite the accompanying paper (citation to be added upon publication).

---

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE) for details.

---

## Contact

Nicolas Largotta — largottn@kean.edu
Levine Lab, Kean University
