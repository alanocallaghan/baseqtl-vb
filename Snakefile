rule all:
    input:
        # "vb_mcmc_estimates_all.png",
        # "optimizing_mcmc_estimates_all.png",
        "rds/GT/vb/all.rds",
        "rds/GT/sampling/all.rds",
        "rds/GT/optimizing/all.rds",
        "rds/noGT/vb/all.rds",
        "rds/noGT/sampling/all.rds",
        "rds/noGT/optimizing/all.rds"

rule run:
    threads: 16
    resources: runtime="10:00:00"
    input: "src/run-{type}.R"
    output: "{type}/{method}/all.rds"
    shell:
        """
        module load R/4.0.3
        Rscript {input} -m {wildcards.method}
        """

rule plot:
    input:
        "{method}/all.rds"
    output:
        "{method}_mcmc_estimates_all.png"
    shell:
        """
        module load R/4.0.3
        Rscript ./plot.R -m {wildcards.method}
        """
