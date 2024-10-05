### page.smk --- Snakemake rules to prepare PAGE data

## Filter variants by R^2 and MAF
rule filter_PAGE:
  output:
    vcf="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/chr{CHROM}.vcf.gz",
    index="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/chr{CHROM}.vcf.gz.tbi",
  container: "docker://mgibio/bcftools:1.12"
  threads: 20
  resources:
    runtime="1d"
  shell:
    """
    bcftools filter -i '(INFO/R2 > 0.6 & INFO/MAF > 0.01)' \
      --threads {threads} \
      /proj/yunligrp/users/quansun/PAGE_subset/MEGA_{wildcards.PAGE_SUB}/chr{wildcards.CHROM}.filtered.vcf.gz \
      -Oz -o {output}
    bcftools index -t {output.vcf}
    """

## Get PAGE phenotypes and write to files
rule get_pheno:
  output:
    expand("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/adjusted/pheno_adj_{PAGE_SUB}.txt", PAGE_SUB=PAGE_SUB),
  resources:
    mem_mb=2000
  container: "docker://frankpo/r_projects:0.0.5"
  shell:
    """
    Rscript workflow/scripts/prepare_phenotype.R
    """

## Construct sample folds for cross-validation
rule make_folds:
  input:
    expand("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/adjusted/pheno_adj_{PAGE_SUB}.txt", PAGE_SUB=PAGE_SUB)
  output:
    expand("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_test.txt", FOLD=FOLD, PHENO=PHENO, PAGE_SUB=PAGE_SUB),
    expand("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_validation.txt", FOLD=FOLD, PHENO=PHENO, PAGE_SUB=PAGE_SUB),
    expand("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_train.txt", FOLD=FOLD, PHENO=PHENO, PAGE_SUB=PAGE_SUB)
  container: "docker://frankpo/r_projects:0.0.5"
  shell:
    """
    Rscript workflow/scripts/make_folds.R
    """


## Combine chromosomes
rule combine_chrom:
  input:
    expand("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{{PAGE_SUB}}/chr{CHROM}.vcf.gz", CHROM=CHROM)
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/all_chrom.vcf.gz"
  threads: 20
  container: "docker://mgibio/bcftools:1.12"
  resources:
    runtime="1d"
  shell:
    """
    bcftools concat {input} --threads {threads} \
      -Oz -o {output}
    """

## Convert vcf to plink format
rule make_plink:
  input:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/all_chrom.vcf.gz"
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/all_chrom.bed",
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/all_chrom.bim",
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/all_chrom.fam"
  container: "docker://quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
  threads: 50
  resources:
    runtime="2d",
    mem_mb=50000
  shell:
    """
    plink2 \
      --threads {threads} \
      --memory 50000 \
      --vcf {input} \
      --double-id \
      --keep-allele-order \
      --allow-extra-chr \
      --make-bed \
      --out /work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{wildcards.PAGE_SUB}/all_chrom
    """
