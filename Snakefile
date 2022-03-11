import json

gt_dict = json.load(open("GT_dict.json"))
nogt_dict = json.load(open("noGT_dict.json"))


rule all:
    input:
        "fig/GT/estimates/pt-95.png",
        "fig/noGT/estimates/pt-95.png"

# rule run:
#     threads: 16
#     resources: runtime="10:00:00"
#     input: "src/run-{type}.R"
#     output: "{type}/{method}/all.rds"
#     shell:
#         """
#         module load R/4.0.3
#         Rscript {input} -m {wildcards.method}
#         """

rule plots:
    input:
        [
            expand(
                "rds/GT/{method}/{gene}_{snp}.rds",
                method = ["sampling"],
                gene = gene,
                snp = gt_dict[gene]
            ) for gene in gt_dict.keys()
        ],
        [
            expand(
                "rds/noGT/{method}/{gene}_{snp}_{condition}_skin.rds",
                method = ["sampling"],
                gene = gene,
                condition = ["Psoriasis", "normal"],
                snp = nogt_dict[gene]
            ) for gene in nogt_dict.keys()
        ]
    output:
        "fig/GT/estimates/pt-95.png",
        "fig/noGT/estimates/pt-95.png"
    shell:
        """
        module load R/4.0.3
        Rscript ./plot-times.R -m GT
        Rscript ./plot-times.R -m noGT
        """


rule run_gt:
    input: "GT_dict.json"
    resources: runtime="01:00:00"
    output: "rds/GT/{method}/{gene}_{snp}.rds"
    shell:
        """
        module load R/4.0.3
        Rscript src/run-snp.R -m {wildcards.method} -g "{wildcards.gene}" -s "{wildcards.snp}"
        """


rule run_nogt:
    input: "noGT_dict.json"
    resources: runtime="01:00:00"
    output:
        "rds/noGT/{method}/{gene}_{snp}_normal_skin.rds",
        "rds/noGT/{method}/{gene}_{snp}_Psoriasis_skin.rds",
    shell:
        """
        module load R/4.0.3
        Rscript src/run-snp.R -n -m {wildcards.method} -g "{wildcards.gene}" -s "{wildcards.snp}"
        """

# rule plot:
#     input:
#         "{method}/all.rds"
#     output:
#         "{method}_mcmc_estimates_all.png"
#     shell:
#         """
#         module load R/4.0.3
#         Rscript ./plot.R -m {wildcards.method}
#         """
