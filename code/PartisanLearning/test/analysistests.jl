import Pkg
Pkg.activate("../")
JULIA_PYTHONCALL_EXE = "/home/marcelovmaciel/miniconda3/bin/python"

using Revise
import PartisanLearning as pl
const pla = pl.PartyLabel
import CSV
# PythonCall.Deps.add(conda_channels = ["conda-forge"],
#                     conda_packages = ["SALib"]
#                     )
