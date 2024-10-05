### ukb.smk -- Snakemake rules to performing clumping on UKB GWAS summary statistics

## Download UKB summary statistics
rule get_sumstats:
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PHENO}.tsv.bgz"
  params:
    url = lambda wildcards: links[wildcards.PHENO]
  container: "docker://mgibio/bcftools:1.12"
  shell:
    """
    ls /work/users/f/r/frocko/
    wget -O {output} {params.url}
    """

## Format UKB summary statistics for PLINK
rule make_plink_assoc:
  input:
    sumstat_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PHENO}.tsv.bgz",
    bim_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/all_chrom.bim"
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PHENO}_{PAGE_SUB}.assoc"
  container: "docker://frankpo/r_projects:0.0.5"
  resources:
    mem_mb = 30000,
    run_time = "2h"
  shell:
    """
    Rscript workflow/scripts/make_plink_assoc.R \
      --input {input.sumstat_file} --output {output} \
      --bim_file {input.bim_file}
    """

## Clump SNPs
rule clump_assoc:
  input:
    assoc="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PHENO}_{PAGE_SUB}.assoc",
    bim_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/all_chrom.bim",
    bed_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/all_chrom.bed",
    fam_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/all_chrom.fam"
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PAGE_SUB}_{PHENO}.clumped"
  container: "docker://quay.io/biocontainers/plink:1.90b6.21--hec16e2b_4"
  resources:
    runtime = "3h",
    mem_mb = 30000
  shell:
    """
    plink \
      --bfile /work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{wildcards.PAGE_SUB}/all_chrom \
      --clump {input.assoc} \
      --clump-p1 0.05 \
      --clump-p2 0.05 \
      --clump-kb 250 \
      --clump-r2 0.1 \
      --out /work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{wildcards.PAGE_SUB}_{wildcards.PHENO}
    """

## List clumped SNPs
rule get_clumped_snps:
  input:
    clump_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PAGE_SUB}_{PHENO}.clumped"
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PAGE_SUB}_{PHENO}_snps.txt"
  container: "docker://mgibio/bcftools:1.12"
  shell:
    """
    awk '{{print $3}}' {input.clump_file} > {output}
    """
