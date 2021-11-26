import Pkg
Pkg.activate("../")
JULIA_PYTHONCALL_EXE = "/home/marcelovmaciel/miniconda3/bin/python"

using Revise
import PartisanLearning as pl
const pla = pl.PartyLabel
import CSV
using Base.Filesystem

# PythonCall.Deps.add(conda_channels = ["conda-forge"],
#                     conda_packages = ["SALib"]
#                     )

results = [pla.get_ParametizationMeasuresMeans(p) for p in 1:1024]

output_matrix = vcat(results...)



CSV.write(joinpath(pla.datapath,"output_matrix.csv"),
          output_matrix)
