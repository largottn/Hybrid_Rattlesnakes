#!/usr/bin/env bash
###############################################################################
# admixture_filtering.sh
#
# Read processing, variant calling, and variant filtering pipeline used to
# prepare whole-genome and ddRAD samples for ADMIXTURE ancestry analysis.
#
# -----------------------------------------------------------------------------
# SOURCE / ATTRIBUTION
#   This pipeline was developed by the Castoe Lab (University of Texas at
#   Arlington) and was provided to me by one of the authors. It corresponds to
#   the variant-calling and filtering workflow described in:
#
#     Maag, D. W., Francioli, Y. Z., Shaw, N., Soni, A. Y., Castoe, T. A.,
#     Schuett, G. W., & Clark, R. W. (2023). Hunting behavior and feeding
#     ecology of Mojave rattlesnakes (Crotalus scutulatus), prairie
#     rattlesnakes (Crotalus viridis), and their hybrids in southwestern
#     New Mexico. Ecology and Evolution, 13, e10683.
#     https://doi.org/10.1002/ece3.10683
#
#   All credit for the pipeline design belongs to the original authors. It is
#   reproduced here with attribution and formatting changes only.
#
# -----------------------------------------------------------------------------
# REFERENCE GENOME
#   Crotalus viridis genome (CroVir_genome_L77pg_16Aug2017.final_rename.fasta;
#   Schield et al. 2019). Must be indexed for BWA (bwa index), samtools
#   (samtools faidx), and GATK (a .dict via Picard CreateSequenceDictionary).
#
# SOFTWARE (versions as reported in Maag et al. 2023 where available)
#   Trimmomatic v0.39    BWA v0.7.17         samtools
#   Picard v2.22.6       GATK v4.1.9.0       VCFtools v0.1.17
#   bcftools / htslib (tabix)                ADMIXTURE v1.3.0 (downstream)
#
# OVERVIEW
#   PART 1 (per sample): trim -> map -> add read groups -> mark duplicates ->
#                        call variants (GATK HaplotypeCaller, GVCF mode)
#   PART 2 (joint):      combine GVCFs -> joint genotype -> hard filter ->
#                        select biallelic SNPs/INDELs -> coverage, missingness,
#                        MAF, and sex-chromosome filters -> final VCF
#   The final VCF (${finalfile}_7.g.vcf.gz) is the input for ADMIXTURE.
#
# EXECUTION NOTE
#   This is a documented, step-by-step pipeline, not a turnkey script. Part 1
#   is run once per sample (e.g. in a loop or array job). Part 2 is run once on
#   the full sample set. One step (97.5th percentile of coverage, step 2e) is
#   computed and then entered manually before step 2h. Fill in all paths and
#   thresholds in the USER-DEFINED VARIABLES block before running.
#
###############################################################################


###############################################################################
# USER-DEFINED VARIABLES  --  edit these before running
###############################################################################

# Reference genome (indexed for BWA, samtools, and GATK; see REFERENCE GENOME)
reference="/path/to/CroVir_genome_L77pg_16Aug2017.final_rename.fasta"

# Sample ID for the per-sample steps in PART 1
# (run PART 1 once per sample, e.g. by looping over a list of IDs)
i="sampleID"

# GATK sample-name map for GenomicsDBImport
#   tab-delimited, two columns: sample_name <TAB> /path/to/${i}.raw.snps.indels.g.vcf.gz
files_to_process="/path/to/mapfile.map"

# Base name for the joint (population) VCF outputs
finalfile="name_of_final_vcf"

# BED / intervals file of repeat elements to mask during VariantFiltration
mask="/path/to/repeat/file.bed"

# File listing female sample IDs (one per line), used for the Z-chromosome filter
females="Female.IDs.txt"

# Minor allele frequency cutoff (Maag et al. 2023 used 0.05)
maf="0.05"

# 97.5th percentile of mean site depth.
# Leave blank for now: it is CALCULATED in step 2e, then pasted back in here
# before running the high-coverage filter in step 2h.
coverage97percentile=""

# Threads
threads=8

# Paths to Java tool jars (adjust to your install; conda installs often provide
# `trimmomatic` and `picard` wrapper commands instead of bare jars)
trimmomatic_jar="/path/to/trimmomatic.jar"
picard_jar="/path/to/picard.jar"


###############################################################################
# PART 1  --  PER-SAMPLE READ PROCESSING AND VARIANT CALLING
# Run once per sample (loop over sample IDs and set $i each time).
#
# Expected directory layout (create these before running):
#   1.RAW.fastqs/                  raw paired FASTQs
#   2.PE.trim.fastqs/              trimmed, paired reads
#   3.PE.trim.unpaired.fastqs/     trimmed, unpaired reads
#   4.RAW.bam.files/               sorted BAMs
#   5.readgroups.added.bam/        BAMs with read groups
#   6.marked.duplicate.bam.files/  duplicate-marked BAMs
#   tmp_bamfiles/  temp_haplotype/ scratch
###############################################################################

## 1a. Trim reads (Trimmomatic, paired-end)
java -Xmx16g -jar "$trimmomatic_jar" PE -phred33 -threads "$threads" \
    1.RAW.fastqs/${i}*R1* 1.RAW.fastqs/${i}*R2* \
    2.PE.trim.fastqs/${i}_R1.P.trim.fq.gz 3.PE.trim.unpaired.fastqs/${i}_R1.U.trim.fq.gz \
    2.PE.trim.fastqs/${i}_R2.P.trim.fq.gz 3.PE.trim.unpaired.fastqs/${i}_R2.U.trim.fq.gz \
    LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30

## 1b. Map trimmed paired reads to the reference (BWA-MEM) and sort (samtools)
bwa mem -t "$threads" "$reference" \
    2.PE.trim.fastqs/${i}_R1.P.trim.fq.gz 2.PE.trim.fastqs/${i}_R2.P.trim.fq.gz \
    | samtools sort -@ "$threads" -O bam -T tmp_bamfiles/temp_${i} -o 4.RAW.bam.files/${i}.bam
samtools index 4.RAW.bam.files/${i}.bam

## 1c. Add read groups (Picard AddOrReplaceReadGroups)
##   (platform, library, etc.) to match your sequencing run.
java -Xmx16g -jar "$picard_jar" AddOrReplaceReadGroups \
    I=4.RAW.bam.files/${i}.bam \
    O=5.readgroups.added.bam/${i}.rg.add.bam \
    RGID=${i} RGLB=lib_${i} RGPL=ILLUMINA RGPU=unit_${i} RGSM=${i}
samtools index 5.readgroups.added.bam/${i}.rg.add.bam

## 1d. Mark duplicate reads (Picard MarkDuplicates; duplicates flagged, not removed)
java -Xmx16g -jar "$picard_jar" MarkDuplicates \
    SORTING_COLLECTION_SIZE_RATIO=0.1 \
    I=5.readgroups.added.bam/${i}.rg.add.bam \
    O=6.marked.duplicate.bam.files/${i}.rg.add.md.bam \
    REMOVE_DUPLICATES=false \
    M=6.marked.duplicate.bam.files/${i}_marked_dup_metrics.txt
samtools index 6.marked.duplicate.bam.files/${i}.rg.add.md.bam

## 1e. Call variants per sample (GATK HaplotypeCaller, GVCF mode)
gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R "$reference" \
    --native-pair-hmm-threads 2 \
    --verbosity ERROR \
    --ERC GVCF \
    --output-mode EMIT_ALL_CONFIDENT_SITES \
    -I 6.marked.duplicate.bam.files/${i}.rg.add.md.bam \
    -O temp_haplotype/${i}.raw.snps.indels.g.vcf.gz


###############################################################################
# PART 2  --  JOINT GENOTYPING AND FILTERING
# Run once on the full set of per-sample GVCFs.
###############################################################################

## 2a. Build a GenomicsDB from all per-sample GVCFs
##   $files_to_process is a tab-delimited map of sample_name -> GVCF path.
gatk --java-options "-Djava.io.tmpdir=tmp_java -Xms4G -Xmx4G -XX:ParallelGCThreads=2" GenomicsDBImport \
    --genomicsdb-workspace-path gdb \
    -R "$reference" \
    --sample-name-map "$files_to_process" \
    --tmp-dir tmp \
    --max-num-intervals-to-import-in-parallel 3

## 2b. Joint genotyping across all samples (GATK GenotypeGVCFs)
gatk --java-options "-Xms8G -Xmx8G -XX:ParallelGCThreads=2" GenotypeGVCFs \
    -R "$reference" \
    -V gendb://gdb \
    -O ${finalfile}_1.g.vcf.gz

## 2c. Hard-filter variants (GATK VariantFiltration) + flag repeat-masked sites
gatk --java-options "-Xmx4g -Xms4g" VariantFiltration \
    --cluster-size 3 \
    --cluster-window-size 10 \
    -V ${finalfile}_1.g.vcf.gz \
    -filter "QD < 2.0"              --filter-name "QD2" \
    -filter "QUAL < 30.0"           --filter-name "QUAL30" \
    -filter "SOR > 3.0"             --filter-name "SOR3" \
    -filter "FS > 60.0"             --filter-name "FS60" \
    -filter "MQ < 40.0"             --filter-name "MQ40" \
    -filter "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -mask "$mask" \
    --mask-name REP \
    --verbosity ERROR \
    -O ${finalfile}_2.g.vcf.gz

## 2d. Keep only biallelic SNPs and INDELs that passed all filters (GATK SelectVariants)
gatk --java-options "-Xmx4g -Xms4g" SelectVariants \
    -R "$reference" \
    -V ${finalfile}_2.g.vcf.gz \
    --select-type-to-include SNP \
    --select-type-to-include INDEL \
    --select-type-to-exclude MIXED \
    --select-type-to-exclude SYMBOLIC \
    --select-type-to-exclude MNP \
    --restrict-alleles-to BIALLELIC \
    --exclude-filtered \
    -O ${finalfile}_3.g.vcf.gz

## 2e. Coverage percentile  --  *** MANUAL STEP ***
##   Compute the 97.5th percentile of mean site depth, then paste the printed
##   value into $coverage97percentile at the top of this script before step 2h.
file=${finalfile}_3.g.vcf.gz
filename=${file##*/}

vcftools --gzvcf "$file" --site-mean-depth --out coverage_$filename

echo "97.5 percentile of coverage::"
sort -k3 -n coverage_$filename.ldepth.mean | awk '{all[NR] = $3} END{print all[int(NR*0.975)]}'

rm coverage_$filename.ldepth.mean
rm coverage_$filename.log

## 2f. Build a list of female-heterozygous Z-chromosome sites to exclude
output="Female.HET.list.kept.sites.exclude.recent.stratum.PAR"

# Pull Z-chromosome genotypes for females
vcftools --gzvcf ${finalfile}_3.g.vcf.gz --chr scaffold-Z --keep "$females" \
    --recode --recode-INFO-all --out TEMP.HET.Z.FEMALE.snps.g

# Keep sites where more than one female is heterozygous
bcftools view --threads "$threads" -i 'COUNT(GT="het") > 1' \
    -O z -o TEMP2.HET.Z.FEMALE.snps.g.vcf.gz TEMP.HET.Z.FEMALE.snps.g.recode.vcf

# Get the list of those sites
vcftools --gzvcf TEMP2.HET.Z.FEMALE.snps.g.vcf.gz --kept-sites --out Female.HET.list

# Drop the recent stratum and PAR (keep positions <= 98,000,000)
awk '{ if ($2 <= 98000000) { print } }' Female.HET.list.kept.sites > "$output"

## 2g. Set zero-coverage genotypes to missing (bcftools +setGT)
bcftools +setGT ${finalfile}_3.g.vcf.gz -- -t q . -n ./. -i 'FMT/DP==0' \
    > ${finalfile}_4.g.vcf.gz
tabix -p vcf ${finalfile}_4.g.vcf.gz

## 2h. Set very high-coverage genotypes to missing (bcftools filter)
##   Uses the 97.5th-percentile value computed in step 2e.
bcftools filter --threads 4 -e "FORMAT/DP > $coverage97percentile" \
    --set-GTs . -O z -o ${finalfile}_5.g.vcf.gz \
    ${finalfile}_4.g.vcf.gz
tabix -p vcf ${finalfile}_5.g.vcf.gz

## 2i. Exclude the female-heterozygous Z-chromosome sites
vcftools --gzvcf ${finalfile}_5.g.vcf.gz \
    --exclude-positions "$output" \
    --recode --recode-INFO-all \
    --stdout | gzip > ${finalfile}_6.g.vcf.gz

## 2j. Final filter: max 20% missing, biallelic, MAF cutoff
vcftools --gzvcf ${finalfile}_6.g.vcf.gz \
    --max-missing 0.8 --min-alleles 2 --max-alleles 2 --maf "$maf" \
    --recode --stdout | gzip > ${finalfile}_7.g.vcf.gz


###############################################################################
# DONE.
###############################################################################
