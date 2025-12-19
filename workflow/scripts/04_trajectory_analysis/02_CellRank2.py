
! source activate /bigdata/godziklab/shared/Xinru/cellrank_env 


import numpy as np

import cellrank as cr
import scanpy as sc
import scvelo as scv

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2

import warnings

warnings.simplefilter("ignore", category=UserWarning)


adata = scv.read("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/velocyto/EOS_scVelo_adata_dynamics.h5ad")
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

scv.pp.moments(adata)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

