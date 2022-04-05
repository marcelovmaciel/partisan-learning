 import Pkg
 Pkg.activate("../")
 JULIA_PYTHONCALL_EXE = "/home/marcelovmaciel/miniconda3/bin/python"

 using Revise
 import PartisanLearning as pl
 const pl.= pl.PartyLabel
import CSV

# PythonCall.Deps.add(conda_channels = ["conda-forge"],
#                     conda_packages = ["SALib"]
#                     )


varnames =  ["ncandidates",  "κ", "δ"]
bounds = [[2.,15.], [0.,7.], [0.5,7.] ]


design_matrix = pl.boundsdict_toparamsdf(varnames, bounds)

design_matrix[!, :ncandidates] = Int.(round.(design_matrix[!, :ncandidates]))

CSV.write("../../../data/saltelli_matrix.csv", design_matrix)
