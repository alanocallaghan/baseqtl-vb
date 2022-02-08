#!/usr/bin/env python3

import theano.tensor as tt
import pymc3 as pm
import pymc3.math as pmath
import numpy
import json

import rpy2.robjects as robjects
import time

from joblib import Parallel, delayed
    
def converter(all):
    da = dict(all.items())
    for k in da.keys():
        da[k] = geneconverter(da[k])
    return da

def geneconverter(snps):
    ds = dict(snps.items())
    for k in ds.keys():
        ds[k] = subconverter(ds[k])
    return ds


def subconverter(snp): 
    d = dict(snp.items())
    return({
        "g": numpy.array(d["g"]), "Y": numpy.array(d["Y"]), "cov": numpy.array(d["cov"])
    })


def gt_nb(x):
    
    Y = x["Y"]
    g = x["g"]
    cov = x["cov"]
    cov = [x[2] for x in cov]
    cov = numpy.array(cov)
    g = numpy.array(g)
    Y = numpy.array(Y)

    with pm.Model() as model:
        N = len(Y)
        mu = numpy.array([0, 0])
        sigma = numpy.array([0.0309, 0.3479])
        weights = numpy.array([0.9736, 0.0264])
        bj = pm.NormalMixture("bj", w = weights, mu = mu, sigma = sigma)
        l1pebj = pmath.log(pmath.exp(bj) + 1) - pmath.log(2)
        phi = pm.Gamma("phi", 1, 0.01)
        alpha = pm.Normal("alpha", 6, 4)
        beta = pm.Normal("beta", 0, 2.5)
        lmu = alpha + (cov * beta)
        for i in range(1, N):
            if abs(g[i] == 1):
                tt.set_subtensor(lmu[i], lmu[i] + l1pebj)
            if g[i] == 2:
                tt.set_subtensor(lmu[i], lmu[i] + bj)
        obs = pm.NegativeBinomial(
            "obs",
            mu = pmath.exp(lmu),
            alpha = phi,
            observed = Y,
            shape = N
        )
        fit = pm.fit(method = "advi")
        # samples = fit.sample(1000)
        # [samples.get_values(x) for x in ('alpha', 'beta', 'bj', 'phi')]
        # samples.get_values("bj")
    return({
            "location": fit.bij.rmap(fit.mean.eval())["bj"].item(),
            "scale": fit.bij.rmap(fit.std.eval())["bj"].item()
    })


readRDS = robjects.r['readRDS']
input = readRDS("rds/gt_data.rds")

input = converter(input)

# out = dict()
# for gene in input.keys():
#     print(gene)
#     gene_data = input[gene]
#     for snp in gene_data.keys():
#         print(snp)
#         snp_data = gene_data[snp]        
#         tic = time.process_time()
#         post = gt_nb(snp_data)
#         toc = time.process_time()
#         out[gene + "_" + snp] = {"samples": post, "time": toc - tic}
#         with open("pymc_samples/" + gene + "_" + snp + ".json", "w") as outfile:
#             json.dump(out[gene + "_" + snp], outfile)

def process(snp):
    snp_data = gene_data[snp]        
    tic = time.process_time()
    post = gt_nb(snp_data)
    toc = time.process_time()
    out = {"samples": post, "time": toc - tic}
    with open("pymc_samples/" + gene + "_" + snp + ".json", "w") as outfile:
        json.dump(out, outfile)


for gene in input.keys():
    print(gene)
    gene_data = input[gene]
    Parallel(n_jobs=16)(delayed(process)(snp) for snp in gene_data.keys())


# with open("pymc_samples.json", "w") as outfile:
#     json.dump(out, outfile)

# gene = "ENSG00000015475"
# snp = "17718119:G:T"
