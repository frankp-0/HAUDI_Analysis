### fbm.smk --- Snakemake rules for constructing the FBMs

rule make_target_regions:
  input:
    flare="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/{PAGE_SUB}/recoded.chr{CHROM}.anc.snps.txt",
    page="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PAGE_SUB}_{PHENO}_snps.txt"
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_regions.txt",
  container: "docker://frankpo/r_projects:0.0.5"
  resources:
    runtime="1h"
  shell:
    """
    Rscript workflow/scripts/make_target_regions.R \
      --flare_snps_file {input.flare} --clumped_snps_file {input.page} --output {output}
    """

rule make_flare_subset_vcf:
  input:
    flare="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/{PAGE_SUB}/chr{CHROM}.anc.vcf.gz",
    regions="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_regions.txt"
  output:
    vcf=temp("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_subset.anc.vcf.gz"),
    index=temp("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_subset.anc.vcf.gz.tbi")
  container: "docker://mgibio/bcftools:1.12"
  resources:
    runtime="2h"
  shell:
    """
    bcftools view -R {input.regions} {input.flare} -Oz -o {output.vcf}
    bcftools index -t {output.vcf}
    """


rule make_page_subset_vcf:
  input:
    page="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/chr{CHROM}.vcf.gz",
    regions="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_regions.txt"
  output:
    vcf=temp("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_subset.vcf.gz"),
    index=temp("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_subset.vcf.gz.tbi")
  container: "docker://mgibio/bcftools:1.12"
  resources:
    runtime="2h"
  shell:
    """
    bcftools view -R {input.regions} {input.page} -Oz -o {output.vcf}
    bcftools index -t {output.vcf}
    """

rule make_lanc_vcf:
  input:
    flare="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_subset.anc.vcf.gz",
    flare_index="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_subset.anc.vcf.gz.tbi",
    page="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_subset.vcf.gz",
    page_index="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_fbm_subset.vcf.gz.tbi"
  output:
    vcf=temp("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_lanc.vcf.gz"),
    index=temp("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/chr{CHROM}_{PHENO}_lanc.vcf.gz.tbi")
  resources:
    runtime="3h"
  container: "docker://mgibio/bcftools:1.12"
  shell:
    """
    bcftools annotate -c FORMAT -a {input.flare} {input.page} -Oz -o {output.vcf}
    bcftools index -t {output.vcf}
    """

rule concat_lanc_vcf:
  input:
    vcf_list=expand("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{{PAGE_SUB}}/chr{CHROM}_{{PHENO}}_lanc.vcf.gz", CHROM=CHROM),
    index_list=expand("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{{PAGE_SUB}}/chr{CHROM}_{{PHENO}}_lanc.vcf.gz.tbi", CHROM=CHROM),
  output:
    vcf="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/all_chrom_{PHENO}_lanc.vcf.gz",
    index="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/all_chrom_{PHENO}_lanc.vcf.gz.tbi",
  resources:
    runtime="3h"
  container: "docker://mgibio/bcftools:1.12"
  shell:
    """
    bcftools concat {input.vcf_list} -Oz -o {output.vcf}
    bcftools index -t {output.vcf}
    """

rule make_fbm:
  input:
    vcf="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/all_chrom_{PHENO}_lanc.vcf.gz",
    index="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/all_chrom_{PHENO}_lanc.vcf.gz.tbi",
  output:
    rds="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc.rds",
    bk="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc.bk",
    info="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc_info.txt"
  resources:
    runtime="1d",
    mem_mb=100000
  container: "docker://frankpo/r_projects:0.0.5"
  shell:
    """
    Rscript workflow/scripts/run_make_fbm.R \
      --vcf_file {input.vcf} \
      --fbm_pref /work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{wildcards.PAGE_SUB}/fbm/all_chrom_{wildcards.PHENO}_lanc
    """
