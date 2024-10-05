### pgs.smk --- Snakemake rules to conduct PGS analysis

## Fit HAUDI models and estimate results
rule run_haudi:
  input:
    rds="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc.rds",
    bk="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc.bk",
    info="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc_info.txt",
    pheno_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/adjusted/pheno_adj_{PAGE_SUB}.txt",
    train_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_train.txt",
    test_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_test.txt",
    validation_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_validation.txt",
    clump_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PAGE_SUB}_{PHENO}.clumped"
  output:
    model="/work/users/f/r/frocko/HAUDI/PAGE_analysis/models/haudi/{PAGE_SUB}/model_{PHENO}_{FOLD}_{n_snps}_{GAMMA}.rds",
    results="/work/users/f/r/frocko/HAUDI/PAGE_analysis/model_results/haudi/{PAGE_SUB}/results_{PHENO}_{FOLD}_{n_snps}_{GAMMA}.txt"
  container: "docker://frankpo/r_projects:0.0.5"
  resources:
    mem_mb=20000,
    runtime="4h"
  shell:
    """
    Rscript workflow/scripts/run_haudi.R \
      --rds {input.rds} \
      --info {input.info} \
      --pheno_file {input.pheno_file} \
      --train_file {input.train_file} \
      --test_file {input.test_file} \
      --validation_file {input.validation_file} \
      --clump_file {input.clump_file} \
      --samples_file source_data/{wildcards.PAGE_SUB}_samples.txt \
      --n_snps {wildcards.n_snps} \
      --gamma {wildcards.GAMMA} \
      --pheno {wildcards.PHENO} \
      --model {output.model} \
      --results {output.results}
    """

## Fit Lasso models and estimate results
rule run_lasso:
  input:
    rds="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc.rds",
    bk="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc.bk",
    info="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc_info.txt",
    pheno_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/adjusted/pheno_adj_{PAGE_SUB}.txt",
    train_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_train.txt",
    test_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_test.txt",
    validation_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_validation.txt",
    clump_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PAGE_SUB}_{PHENO}.clumped"
  output:
    model="/work/users/f/r/frocko/HAUDI/PAGE_analysis/models/lasso/{PAGE_SUB}/model_{PHENO}_{FOLD}_{n_snps}.rds",
    results="/work/users/f/r/frocko/HAUDI/PAGE_analysis/model_results/lasso/{PAGE_SUB}/results_{PHENO}_{FOLD}_{n_snps}.txt"
  container: "docker://frankpo/r_projects:0.0.5"
  resources:
    mem_mb=15000,
    runtime="4h"
  shell:
    """
    Rscript workflow/scripts/run_lasso.R \
      --rds {input.rds} \
      --info {input.info} \
      --pheno_file {input.pheno_file} \
      --train_file {input.train_file} \
      --test_file {input.test_file} \
      --validation_file {input.validation_file} \
      --clump_file {input.clump_file} \
      --samples_file source_data/{wildcards.PAGE_SUB}_samples.txt \
      --n_snps {wildcards.n_snps} \
      --pheno {wildcards.PHENO} \
      --model {output.model} \
      --results {output.results}
    """

## Fit GAUDI models and estimate results
rule run_gaudi:
  input:
    rds="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc.rds",
    bk="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc.bk",
    info="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/fbm/{PAGE_SUB}/fbm/all_chrom_{PHENO}_lanc_info.txt",
    pheno_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/adjusted/pheno_adj_{PAGE_SUB}.txt",
    train_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_train.txt",
    test_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_test.txt",
    validation_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold{FOLD}_{PHENO}_{PAGE_SUB}_validation.txt",
    clump_file="/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/UKB/{PAGE_SUB}_{PHENO}.clumped"
  output:
    model="/work/users/f/r/frocko/HAUDI/PAGE_analysis/models/gaudi/{PAGE_SUB}/model_{PHENO}_{FOLD}_{n_snps}_{GAMMA}.rds",
    results="/work/users/f/r/frocko/HAUDI/PAGE_analysis/model_results/gaudi/{PAGE_SUB}/results_{PHENO}_{FOLD}_{n_snps}_{GAMMA}.txt"
  container: "docker://frankpo/r_projects:0.0.5"
  resources:
    mem_mb=25000,
    runtime="1d"
  shell:
    """
    Rscript workflow/scripts/run_gaudi.R \
      --rds {input.rds} \
      --info {input.info} \
      --pheno_file {input.pheno_file} \
      --train_file {input.train_file} \
      --test_file {input.test_file} \
      --validation_file {input.validation_file} \
      --clump_file {input.clump_file} \
      --samples_file source_data/{wildcards.PAGE_SUB}_samples.txt \
      --n_snps {wildcards.n_snps} \
      --gamma {wildcards.GAMMA} \
      --pheno {wildcards.PHENO} \
      --model {output.model} \
      --results {output.results}
    """
