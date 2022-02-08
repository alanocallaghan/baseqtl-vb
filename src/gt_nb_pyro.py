import pyro
from pyro.distributions \
    import Gamma, NegativeBinomial, Normal, Cauchy
from pyro import plate, sample, param
from pyro import distributions as dist
from pyro.infer import EmpiricalMarginal, SVI, Trace_ELBO
import torch
from torch.distributions import constraints
from torch import ones, float64, exp, tensor, matmul, diag, zeros, log, tensor, Size
from pyro.contrib.autoguide import AutoDiagonalNormal

from statistics import median
import math
import json

import numpy
import matplotlib.pyplot as plt


import rpy2.robjects as robjects
import time


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

# https://forum.pyro.ai/t/mixture-of-poissons-with-large-sample-size/1901
class MixtureNormal(dist.TorchDistribution):
    support = dist.Normal.support
    has_rsample = False

    def __init__(self, mu, sigma, theta):
        self._normal = dist.Normal(mu, sigma)
        self.theta = theta
        self._num_component = theta.shape[-1]

        cdbs = self._normal.batch_shape[:-1]
        event_shape = self._normal.event_shape
        self._event_ndims = len(event_shape)

        super(MixtureNormal, self).__init__(batch_shape=cdbs,
                                            event_shape=event_shape,
                                            validate_args=False)

    def sample(self, sample_shape=Size()):
        return ones(sample_shape + self.shape())

    def log_prob(self, value):
        contributions = self._normal.log_prob(value) + self.theta.log()

        return contributions.logsumexp(-1)


def gt_nb(Y, cov, g):

    bj = sample("bj", MixtureNormal(
            tensor([0, 0]),
            tensor([0.0309, 0.3479]),
            tensor([0.9736, 0.0264])
        )
    )
    phi = sample("phi", Gamma(1, 0.01))
    alpha = sample("alpha", Normal(6, 4))
    beta = sample("beta", Cauchy(6, 4))
    lmu = alpha + (cov * beta)
    l1pebj = math.log(math.exp(bj) + 1) - math.log(2)
    
    for i in range(len(Y)):
        if abs(g[i]) == 1:
            lmu[i] = lmu[i] + l1pebj
        if g[i] == 2:
            lmu[i] = lmu[i] + bj

    r = phi
    p = exp(lmu) / (exp(lmu) + phi)
    sample(
        "counts", 
        NegativeBinomial(total_count = r, probs = p).to_event(1), 
        obs = Y
    )


def posterior_samples_svi(svi, sites, **kwargs):
    posterior_samples = {
        site: EmpiricalMarginal(posterior, sites = site)\
            .sample((1000, )).cpu().detach().numpy()
        for site in sites}
    return(posterior_samples)


def run_svi(
        svi, 
        num_iters = 100000, 
        window = 100, 
        threshold = 1e-5,
        **kwargs
        ):

    elbo_h = None
    norms = numpy.ones(num_iters)
    elbos = numpy.zeros(num_iters)
    for step in range(num_iters):
        elbo = svi.step(**kwargs)
        elbos[step] = elbo
        if not elbo_h == None:
            norms[step] = abs((elbo - elbo_h) / elbo_h)
            r = range(step - window, step)
            r = [i for i in r if i >= 0]
            mean_norm = sum(norms[r]) / len(r)
            median_norm = median(norms[r])
            if mean_norm < threshold or median_norm < threshold:
                norms = norms[norms != 1]
                elbos = elbos[elbos != 0]
                break
        if (step % 100 == 0):
            print("Iteration: %05d; ELBO: %10.4f; Norm: %04.4f" % (step, elbo, norms[step]))
        elbo_h = elbo
    return svi, elbos, norms


def run_gt_nb(x):

    Y = x["Y"]
    g = x["g"]
    cov = x["cov"]
    cov = [x[2] for x in cov]
    cov = numpy.array(cov)
    g = numpy.array(g)
    Y = numpy.array(Y)
    cov = torch.from_numpy(cov)
    g = torch.from_numpy(g)
    Y = torch.from_numpy(Y)

    # http://pyro.ai/examples/svi_part_iv.html
    initial_lr = 0.001
    gamma = 0.1  # final learning rate will be gamma * initial_lr
    lrd = gamma ** (1 / 10000)
    optim = pyro.optim.ClippedAdam({'lr': initial_lr, 'lrd': lrd})

    guide_gt_nb = AutoDiagonalNormal(gt_nb)
    svi = SVI(
        gt_nb, 
        guide_gt_nb,
        optim = optim,
        loss = Trace_ELBO()
    )
    svi, elbos, norms = run_svi(
        svi,
        Y = Y,
        g = g,
        cov = cov,
        threshold = 1e-2,
        window = 10,
        num_iters = 10000
    )

    # guide_gt_nb.median(Y, g, cov)
    ls = guide_gt_nb._loc_scale()
    location = ls[0][0].item()
    scale = ls[1][0].item()

    return {"location": location, "scale": scale}


readRDS = robjects.r['readRDS']
input = readRDS("rds/gt_data.rds")

input = converter(input)

out = dict()
for gene in input.keys():
    gene_data = input[gene]
    for snp in gene_data.keys():
        print(gene)
        print(snp)
        snp_data = gene_data[snp]        
        while True:
            try:
                tic = time.process_time()
                post = run_gt_nb(snp_data)
                toc = time.process_time()
            except ValueError:
                continue
            break
        out[gene + "_" + snp] = {"params": post, "time": toc - tic}

with open("pyro_params.json", "w") as outfile:
    json.dump(out, outfile)

# gene="ENSG00000015475"
# snp="18219727:A:G"
