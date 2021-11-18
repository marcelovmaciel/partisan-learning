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

designmatrix = CSV.read("../../../data/saltelli_matrix.csv",
                        pla.DF.DataFrame)

simulation_constants = Dict(:nagents => 500, :nissues => 2, :niterations => 1000)


x = eachrow(designmatrix)[1]

pla.run_analysis_onRow(simulation_constants, x, 1, 1, 1 )
