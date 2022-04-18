import json

gt_dict = json.load(open("GT_dict.json"))
nogt_dict = json.load(open("noGT_dict.json"))
components = ["both", "inter", "intra"]
tols = ["1e-02"]
# tols = ["1e-02", "1e-03"]

shell.prefix("module load R/4.0.3; ")

rule all:
    input:
        # "fig/GT/components/estimates/point-estimates.png",
        "fig/GT/estimates/point-estimates-95.png",
        "fig/GT/estimates/sampling_mcmc_all_categorical_95.png",
        "fig/GT/estimates/vb_mcmc_all_categorical_95.png",
        "fig/noGT/estimates/point-estimates-95.png",
        "fig/noGT/estimates/sampling_mcmc_all_categorical_95.png",
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
        "rds/GT/sampling_combined.rds",
        "rds/GT/vb_1e-02_combined.rds",
        # "rds/GT/vb_1e-03_combined.rds",
        "rds/tmp2.rds",
        # "rds/GT/components/sampling_combined.rds",
        # "rds/GT/components/vb_1e-02_combined.rds",
        # "rds/GT/components/vb_1e-03_combined.rds",
        # "rds/tmp.rds" ## just so the rule counts as unfinished
        # ,
        "rds/noGT/sampling_combined.rds",
        "rds/noGT/vb_1e-02_combined.rds",
        "rds/noGT/vb_1e-03_combined.rds"
    output:
        # "fig/GT/components/estimates/point-estimates.png",
        "fig/GT/estimates/point-estimates-95.png",
        "fig/GT/estimates/sampling_mcmc_all_categorical_95.png",
        "fig/GT/estimates/vb_mcmc_all_categorical_95.png",
        "fig/noGT/estimates/point-estimates-95.png",
        "fig/noGT/estimates/sampling_mcmc_all_categorical_95.png",
        "fig/noGT/estimates/vb_mcmc_all_categorical_95.png"
    shell:
        """
        Rscript ./src/plot-times.R -m GT -t 1e-02
        # Rscript ./src/plot-times.R -m GT -t 1e-03
        Rscript ./src/plot.R -m GT -i vb -t 1e-02
        # Rscript ./src/plot.R -m GT -i vb -t 1e-03
        Rscript ./src/plot.R -m GT -i sampling
        # Rscript ./src/plot-components.R -m GT -t 1e-02
        # Rscript ./src/plot-components.R -m GT -t 1e-03

        Rscript ./src/plot-times.R -m noGT -t 1e-02
        # Rscript ./src/plot-times.R -m noGT -t 1e-03
        Rscript ./src/plot.R -m noGT -i vb -t 1e-02
        # Rscript ./src/plot.R -m noGT -i vb -t 1e-03
        Rscript ./src/plot.R -m noGT -i sampling
        """

rule process:
    resources: runtime="02:00:00"
    input:
        expand(
            "rds/GT/vb_{tol}/{gene}_done",
            tol = tols,
            gene = gt_dict.keys()
        ),
        expand(
            "rds/GT/sampling/{gene}_done",
            gene = gt_dict.keys()
        ),
        expand(
            "rds/noGT/sampling/{gene}_done",
            gene = nogt_dict.keys()
        ),
        expand(
            "rds/noGT/vb_{tol}/{gene}_done",
            tol = tols,
            gene = nogt_dict.keys()
        )
    output:
        "rds/GT/sampling_combined.rds",
        "rds/GT/vb_1e-02_combined.rds",
        "rds/GT/vb_1e-03_combined.rds",
        "rds/tmp2.rds",
        "rds/noGT/sampling_combined.rds",
        "rds/noGT/vb_1e-02_combined.rds",
        "rds/noGT/vb_1e-03_combined.rds"
    shell:
        """
        Rscript src/process-files.R -i sampling -m GT
        Rscript src/process-files.R -i vb -m GT -t 1e-02
        # Rscript src/process-files.R -i vb -m GT -t 1e-03
        Rscript src/process-files.R -i sampling -m noGT
        Rscript src/process-files.R -i vb -m noGT -t 1e-02
        # Rscript src/process-files.R -i vb -m noGT -t 1e-03
        """

rule process_components:
    resources: runtime="02:00:00"
    input:
        "rds/GT/components/sampling/{gene}_done.rds",
        "rds/GT/components/vb_{tol}/{gene}_done.rds"
        # [
        #     expand(
        #         "rds/GT/components/vb_{tol}/{gene}_{snp}.rds",
        #         tol = tols,
        #         gene = gene,
        #         snp = gt_dict[gene]
        #     ) for gene in gt_dict.keys()
        # ],
        # [
        #     expand(
        #         "rds/GT/components/sampling/{gene}_{snp}.rds",
        #         gene = gene,
        #         snp = gt_dict[gene]
        #     ) for gene in gt_dict.keys()
        # ]
    output:
        "rds/GT/components/sampling_combined.rds",
        "rds/GT/components/vb_1e-02_combined.rds",
        "rds/GT/components/vb_1e-03_combined.rds"
        # ,
        # "rds/tmp.rds" ## just so the rule counts as unfinished
    shell:
        """
        Rscript src/process-component-files.R -i vb -m GT -t 1e-02
        # Rscript src/process-component-files.R -i vb -m GT -t 1e-03
        Rscript src/process-component-files.R -i sampling -m GT
        """


rule run_gt_sc:
    input: "dicts/{gene}.json"
    threads: 4
    resources: runtime="05:00:00"
    output:
        "rds/GT/components/sampling/{gene}_done"
        # "rds/GT/components/sampling/{gene}_{snp}.rds"
    shell:
        """
        Rscript src/run-snp-sep.R \
            -m GT \
            -i sampling \
            -g "{wildcards.gene}"
        """

rule run_gt_vc:
    resources: runtime="05:00:00"
    output:
        "rds/GT/components/vb_{tol}/{gene}_done"
    shell:
        """
        Rscript src/run-snp-sep.R \
            -m GT \
            -i vb \
            -g "{wildcards.gene}" \
            -t {wildcards.tol}
        """


rule run_gt_s:
    resources: runtime="05:00:00"
    threads: 4
    output:
        "rds/GT/sampling/{gene}_done"
    shell:
        """
        Rscript src/run-snp.R -m GT -i sampling -g "{wildcards.gene}"
        """

rule run_gt_v:
    resources: runtime="05:00:00"
    output:
        # "rds/GT/vb_{tol}/{gene}_{snp}.rds"
        "rds/GT/vb_{tol}/{gene}_done"
    shell:
        """
        Rscript src/run-snp.R \
            -m GT \
            -i vb \
            -g "{wildcards.gene}" \
            -t {wildcards.tol}
        """

rule run_nogt_s:
    resources: runtime="05:00:00"
    threads: 4
    output:
        "rds/noGT/sampling/{gene}_done"
        # "rds/noGT/sampling/{gene}_{snp}_normal_skin.rds",
        # "rds/noGT/sampling/{gene}_{snp}_Psoriasis_skin.rds"
    shell:
        """
        Rscript src/run-snp.R \
            -m noGT \
            -i sampling \
            -g "{wildcards.gene}"
        """ 

rule run_nogt_v:
    resources: runtime="05:00:00"
    output:
        "rds/noGT/vb_{tol}/{gene}_done"
        # "rds/noGT/vb_{tol}/{gene}_{snp}_normal_skin.rds",
        # "rds/noGT/vb_{tol}/{gene}_{snp}_Psoriasis_skin.rds"
    shell:
        """
        Rscript src/run-snp.R \
            -m noGT \
            -i vb \
            -g "{wildcards.gene}" \
            -t {wildcards.tol}
        """


rule run_lm:
    resources: runtime="01:00:00"
    output:
        "rds/GT/lm-filtering.rds"
    shell:
        """
        Rscript src/run-lm.R
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
