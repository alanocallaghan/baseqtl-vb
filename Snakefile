import json

gt_dict = json.load(open("GT_dict.json"))
nogt_dict = json.load(open("noGT_dict.json"))

shell.prefix("module load R/4.0.3; ")

rule all:
    input:
        "fig/GT/estimates/point-estimates-95.png",
        "fig/noGT/estimates/point-estimates-95.png",
        "fig/GT/estimates/sampling_mcmc_all_categorical_95.png",
        "fig/noGT/estimates/sampling_mcmc_all_categorical_95.png",
        "fig/GT/estimates/vb_mcmc_all_categorical_95.png",
        "fig/noGT/estimates/vb_mcmc_all_categorical_95.png"

# rule run:
#     threads: 16
#     resources: runtime="10:00:00"
#     input: "src/run-{type}.R"
#     output: "{type}/{method}/all.rds"
#     shell:
#         """
#         Rscript {input} -m {wildcards.method}
#         """

rule plots:
    resources: runtime="02:00:00"
    input:
        # [
        #     expand(
        #         "rds/GT/{method}/{gene}_{snp}.rds",
        #         method = ["vb", "sampling"],
        #         gene = gene,
        #         snp = gt_dict[gene]
        #     ) for gene in gt_dict.keys()
        # ],
        # [
        #     expand(
        #         "rds/noGT/{method}/{gene}_{snp}_{condition}_skin.rds",
        #         method = ["vb", "sampling"],
        #         gene = gene,
        #         condition = ["Psoriasis", "normal"],
        #         snp = nogt_dict[gene]
        #     ) for gene in nogt_dict.keys()
        # ]
    output:
        "fig/GT/estimates/point-estimates-95.png",
        "fig/noGT/estimates/point-estimates-95.png",
        "fig/GT/estimates/sampling_mcmc_all_categorical_95.png",
        "fig/noGT/estimates/sampling_mcmc_all_categorical_95.png",
        "fig/GT/estimates/vb_mcmc_all_categorical_95.png",
        "fig/noGT/estimates/vb_mcmc_all_categorical_95.png"
    shell:
        """
        Rscript ./src/plot-times.R -m GT
        Rscript ./src/plot-times.R -m noGT
        Rscript ./src/plot.R -m GT -i vb
        Rscript ./src/plot.R -m GT -i sampling
        Rscript ./src/plot.R -m noGT -i sampling
        Rscript ./src/plot.R -m noGT -i vb
        """

rule process:
    resources: runtime="02:00:00"
    input:
        # [
        #     expand(
        #         "rds/GT/{method}/{gene}_{snp}.rds",
        #         method = ["vb", "sampling"],
        #         gene = gene,
        #         snp = gt_dict[gene]
        #     ) for gene in gt_dict.keys()
        # ],
        # [
        #     expand(
        #         "rds/noGT/{method}/{gene}_{snp}_{condition}_skin.rds",
        #         method = ["vb", "sampling"],
        #         gene = gene,
        #         condition = ["Psoriasis", "normal"],
        #         snp = nogt_dict[gene]
        #     ) for gene in nogt_dict.keys()
        # ]
    output:
        "rds/GT/sampling_combined.rds",
        "rds/GT/vb_combined.rds",
        "rds/noGT/sampling_combined.rds",
        "rds/noGT/vb_combined.rds"
    shell:
        """
        Rscript src/process-files.R -i vb -m GT
        Rscript src/process-files.R -i vb -m noGT
        Rscript src/process-files.R -i sampling -m GT
        Rscript src/process-files.R -i sampling -m noGT
        """


rule run_gt:
    input: "GT_dict.json"
    resources: runtime="01:00:00"
    output:
        "rds/GT/{method}/{gene}_{snp}.rds"
    shell:
        """
        Rscript src/run-snp.R -m GT -i {wildcards.method} -g "{wildcards.gene}" -s "{wildcards.snp}"
        """


rule run_nogt:
    input: "noGT_dict.json"
    resources: runtime="01:00:00"
    output:
        "rds/noGT/{method}/{gene}_{snp}_normal_skin.rds",
        "rds/noGT/{method}/{gene}_{snp}_Psoriasis_skin.rds"
    shell:
        """
        Rscript src/run-snp.R -m noGT -i {wildcards.method} -g "{wildcards.gene}" -s "{wildcards.snp}"
        """

# rule plot:
#     input:
#         "{method}/all.rds"
#     output:
#         "{method}_mcmc_estimates_all.png"
#     shell:
#         """
#         Rscript ./plot.R -m {wildcards.method}
#         """

rule dicts:
    input: "src/make-dicts.R"
    output:
        "GT_dict.json",
        "noGT_dict.json"
    shell:
        """
        Rscript {input}
        """
