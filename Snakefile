import json

gt_dict = json.load(open("GT_dict.json"))
nogt_dict = json.load(open("noGT_dict.json"))
dicts = {"GT": gt_dict, "noGT": nogt_dict}
components = ["both", "inter", "intra"]
# tols = ["1e-02"]
tols = ["1e-02", "1e-03"]
seeds = [7, 14, 21, 28, 35, 42]
gt_in_dir = "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT/"
nogt_in_dir = "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"

# methods = ["sampling", "pathfinder", "pathfinder_parallel"]
methods = ["sampling"]
models = ["GT", "noGT"]

shell.prefix(". /etc/profile.d/modules.sh; module load gcc/9; module load R/4.2.0-icelake; ")

wildcard_constraints:
    gene="[^_/]+",
    seed="\d+",
    tol="[^_/]+"

rule all:
    input:
        expand(
            "fig_{tol}/{model}/estimates/point-estimates-both-99.pdf",
            tol = tols,
            model = models
        ),
        expand(
            "fig_{tol}/GT/diag/KS-meanfield-fullrank.pdf",
            tol = tols
        ),
        expand(
            "fig_{tol}/{model}/diag/full-posteriors-slab.pdf",
            tol = tols,
            model = models
        ),
        expand(
            "fig_{tol}/{model}/variability/time-roc-pep.pdf",
            tol = tols,
            model = models
        ),
        expand(
            "fig_{tol}/GT/roc/time_roc_all_99.pdf",
            tol = tols
        )


rule plots:
    resources: runtime="02:00:00", mem_mb=20000
    input:
        "rds/GT/sampling_combined.rds",
        "rds/GT/vb_{tol}_combined.rds",
        "rds/noGT/sampling_combined.rds",
        "rds/noGT/vb_{tol}_combined.rds"
    output:
        "fig_{tol}/GT/estimates/point-estimates-both-99.pdf",
        "fig_{tol}/noGT/estimates/point-estimates-both-99.pdf",
        "fig_{tol}/GT/variability/time-roc-pep.pdf",
        "fig_{tol}/noGT/variability/time-roc-pep.pdf"
    shell:
        """
        Rscript ./src/plot-times.R -m GT -t {wildcards.tol}
        Rscript ./src/plot-times.R -m noGT -t {wildcards.tol}
        Rscript ./src/plot-variability.R -m GT -t {wildcards.tol}
        Rscript ./src/plot-variability.R -m noGT -t {wildcards.tol}
        # Rscript ./src/plot-components.R -m GT -t {wildcards.tol}
        """


rule rerun_discrepant:
    resources: runtime="04:00:00", mem_mb=20000
    input:
        "rds/{model}_discrepancies_vb_{tol}.rds"
    output:
        "rds/{model}_vb_{tol}_rerun_discrepant_full.rds"
    shell:
        """
        Rscript src/rerun-most-discrepant.R -m {wildcards.model} -t {wildcards.tol}
        """


rule plot_discrepant:
    resources: runtime="04:00:00", mem_mb=20000
    input:
        "rds/{model}_vb_{tol}_rerun_discrepant_full.rds"
    output:
        "fig_{tol}/{model}/diag/full-posteriors-slab.pdf",
        "fig_{tol}/{model}/diag/most_discrepant_method_comp_prior.pdf"
    shell:
        """
        Rscript src/plot-discrepancy-posteriors.R -m {wildcards.model} -t {wildcards.tol}
        Rscript src/plot-discrepancy-sources.R -m {wildcards.model} -t {wildcards.tol}
        Rscript src/plot-discrepancy-misc.R -m {wildcards.model} -t {wildcards.tol}
        """


rule plot_meanfield:
    resources: runtime="04:00:00", mem_mb=50000
    threads: 8
    input:
        expand(
            "rds/GT/vb_{{tol}}/meanfield/{gene}_s{seed}.rds",
            gene = dicts["GT"].keys(),
            seed = 42
        )
    output:
        "fig_{tol}/GT/diag/KS-meanfield-fullrank.pdf"
    shell:
        """
        Rscript ./src/plot-meanfield.R -m GT -t {wildcards.tol}
        """


rule process:
    resources: runtime="02:00:00", mem_mb=20000
    threads: 8
    input:
        [
            expand(
                "rds/{model}/vb_{{tol}}/{gene}_s{seed}.rds",
                model = mod,
                gene = dicts[mod].keys(),
                seed = seeds
            )
            for mod in models
        ]
    output:
        "rds/{model}/vb_{tol}_combined.rds"
    shell:
        """
        Rscript src/process-files.R -i vb -m {wildcards.model} -t {wildcards.tol}
        """


rule process_sampling:
    resources: runtime="02:00:00", mem_mb=20000
    threads: 8
    input:
        [
            expand(
                "rds/{model}/sampling/{gene}_s{seed}.rds",
                model = mod,
                gene = dicts[mod].keys(),
                seed = 42
            ) for mod in models
        ]
    output:
        "rds/{model}/sampling_combined.rds"
    shell:
        """
        Rscript src/process-files.R -i sampling -m {wildcards.model}
        """


rule process_components_sampling:
    resources: runtime="02:00:00"
    input:
        "rds/GT/components/sampling/{gene}_s{seed}.rds",
    output:
        "rds/GT/components/sampling_combined.rds"
    shell:
        """
        Rscript src/process-component-files.R -i sampling -m GT -t {wildcards.tol}
        """


rule process_components:
    resources: runtime="02:00:00"
    input:
        "rds/GT/components/vb_{tol}/{gene}_s{seed}.rds"
    output:
        "rds/GT/components/vb_{tol}_combined.rds"
    shell:
        """
        Rscript src/process-component-files.R -i vb -m GT -t {wildcards.tol}
        """


rule run_gt_components:
    resources: runtime="05:00:00", mem_mb=10000
    threads: lambda wildcards: 4 if wildcards.method == "sampling" else 8
    input:
        gt_in_dir + "/rbias.{gene}.GT.stan1.input.rds"
    output:
        "rds/GT/components/{method}_{tol}/{gene}_s{seed}.rds"
    shell:
        """
        Rscript src/run-snp-sep.R \
            -m GT \
            -i {wildcards.method} \
            -g {wildcards.gene} \
            -t {wildcards.tol} \
            -c {threads} \
            -s {wildcards.seed}
        """


rule run_gt:
    resources: runtime="05:00:00", mem_mb=10000
    threads: 8
    input:
        gt_in_dir + "rbias.{gene}.GT.stan1.input.rds"
    output:
        "rds/GT/vb_{tol}/{gene}_s{seed}.rds"
    shell:
        """
        Rscript src/run-snp.R \
            -m GT \
            -i vb \
            -g {wildcards.gene} \
            -t {wildcards.tol} \
            -c {threads} \
            -s {wildcards.seed}
        """


rule run_gt_sampling:
    resources: runtime="05:00:00", mem_mb=10000
    threads: 8
    input:
        gt_in_dir + "rbias.{gene}.GT.stan1.input.rds"
    output:
        "rds/GT/sampling/{gene}_s{seed}.rds"
    shell:
        """
        Rscript src/run-snp.R \
            -m GT \
            -i sampling \
            -g {wildcards.gene} \
            -c {threads} \
            -s {wildcards.seed}
        """


rule run_nogt:
    resources: runtime="05:00:00", mem_mb=20000
    threads: 8
    input:
        normal = nogt_in_dir + "refbias.{gene}.normal_skin.noGT.stan.input.rds",
        Psoriasis = nogt_in_dir + "refbias.{gene}.Psoriasis_skin.noGT.stan.input.rds"
    output:
        "rds/noGT/vb_{tol}/{gene}_s{seed}.rds"
    shell: 
        """
        Rscript src/run-snp.R \
            -m noGT \
            -i vb \
            -g {wildcards.gene} \
            -t {wildcards.tol} \
            -c {threads} \
            -s {wildcards.seed} 
        """


rule run_nogt_sampling:
    resources: runtime="05:00:00", mem_mb=20000
    threads: 8
    input:
        normal = nogt_in_dir + "refbias.{gene}.normal_skin.noGT.stan.input.rds",
        Psoriasis = nogt_in_dir + "refbias.{gene}.Psoriasis_skin.noGT.stan.input.rds"
    output:
        "rds/noGT/sampling/{gene}_s{seed}.rds"
    shell: 
        """
        Rscript src/run-snp.R \
            -m noGT \
            -i sampling \
            -g {wildcards.gene} \
            -c {threads} \
            -s {wildcards.seed} 
        """


rule run_gt_mf:
    resources: runtime="05:00:00", mem_mb=10000
    threads: 8
    input:
        gt_in_dir + "rbias.{gene}.GT.stan1.input.rds"
    output:
        "rds/GT/vb_{tol}/meanfield/{gene}_s{seed}.rds"
    shell:
        """
        Rscript src/run-meanfield-fullrank.R \
            -m GT \
            -g "{wildcards.gene}" \
            -t {wildcards.tol} \
            -c {threads} \
            -s {wildcards.seed}
        """


rule run_nogt_mf:
    resources: runtime="05:00:00", mem_mb=10000
    threads: 8
    input:
        normal = nogt_in_dir + "refbias.{gene}.normal_skin.noGT.stan.input.rds",
        Psoriasis = nogt_in_dir + "refbias.{gene}.Psoriasis_skin.noGT.stan.input.rds"
    output:
        "rds/noGT/{method}_{tol}/meanfield/{gene}_s{seed}.rds"
    shell:
        """
        Rscript src/run-meanfield-fullrank.R \
            -m noGT \
            -g "{wildcards.gene}" \
            -t {wildcards.tol} \
            -c {threads} \
            -s {wildcards.seed}
        """


rule plot_lm:
    resources: runtime="01:00:00"
    threads: 1
    input:
        "rds/GT/lm-filtering.rds"
    output:
        "fig_{tol}/GT/roc/time_roc_all_99.pdf"
    shell:
        """
        Rscript src/plot-lm.R -t {wildcards.tol}
        """


rule run_lm:
    resources: runtime="01:00:00"
    threads: 8
    output:
        "rds/GT/lm-filtering.rds"
    shell:
        """
        Rscript src/run-lm.R
        """
