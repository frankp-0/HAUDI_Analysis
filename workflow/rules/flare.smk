### flare.smk --- Snakemake rules for conducting flare

## Filter by MAF and sort reference panel VCFs
rule filter_sort_index_reference:
  output:
    vcf="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/reference_panel/chr{CHROM}.vcf.gz",
    index="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/reference_panel/chr{CHROM}.vcf.gz.tbi"
  container: "docker://mgibio/bcftools:1.12"
  resources:
    runtime="1h"
  shell:
    """
    bcftools view -t chr{wildcards.CHROM} -q 0.05:minor \
      source_data/reference_panel/reference1.phased.hg38.chr{wildcards.CHROM}.vcf.gz -Ou | \
      bcftools sort \
      -Oz -o {output}
    bcftools index -t {output}
    """

## Run flare
rule flare:
  input:
    target="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/geno/{PAGE_SUB}/chr{CHROM}.vcf.gz",
    ref="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/reference_panel/chr{CHROM}.vcf.gz"
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/{PAGE_SUB}/chr{CHROM}.anc.vcf.gz"
  container: "docker://eclipse-temurin:8u382-b05-jre-focal"
  threads: 30
  resources:
    mem_mb=50000,
    runtime="1d"
  shell:
    """
    java -Xmx50G -jar lib/flare.jar \
      ref={input.ref} \
      ref-panel=source_data/{wildcards.PAGE_SUB}_map.txt \
      min-maf=0 \
      array=true \
      gt={input.target} \
      map=source_data/plink.GRCh38.map \
      nthreads={threads} \
      seed=34991002 \
      out={output}
        """

rule flare_recode_id:
  input:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/{PAGE_SUB}/chr{CHROM}.anc.vcf.gz"
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/{PAGE_SUB}/recoded.chr{CHROM}.anc.vcf.gz"
  container: "docker://mgibio/bcftools:1.12"
  shell:
    """
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' {input} -Oz -o {output}
    """

## Index flare output file        
rule index_flare:
  input:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/{PAGE_SUB}/recoded.chr{CHROM}.anc.vcf.gz"
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/{PAGE_SUB}/recoded.chr{CHROM}.anc.vcf.gz.tbi"
  container: "docker://mgibio/bcftools:1.12"        
  resources:
    runtime="2h"
  shell:
    """
    bcftools index -t {input}
    """

## List SNPs in flare output file
rule flare_snps:
  input:
    vcf="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/{PAGE_SUB}/recoded.chr{CHROM}.anc.vcf.gz",
    index="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/{PAGE_SUB}/recoded.chr{CHROM}.anc.vcf.gz.tbi"
  output:
    "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/flare/{PAGE_SUB}/recoded.chr{CHROM}.anc.snps.txt"
  container: "docker://mgibio/bcftools:1.12"
  shell:
    """
    bcftools query -f '%ID' {input.vcf} > {output}
    """
