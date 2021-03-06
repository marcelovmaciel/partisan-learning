import Pkg
Pkg.activate("../")
JULIA_PYTHONCALL_EXE = "/home/marcelovmaciel/miniconda3/bin/python"

using Revise
import PartisanLearning as pl
const pl.= pl.PartyLabel
import CSV
using Base.Filesystem

# PythonCall.Deps.add(conda_channels = ["conda-forge"],
#                     conda_packages = ["SALib"]
#                     )

results = [pl.get_ParametizationMeasuresMeans(p) for p in 1:2048]

output_matrix = vcat(results...)

output_matrix

CSV.write(joinpath(pl.datapath,"output_matrix.csv"),
          output_matrix)
