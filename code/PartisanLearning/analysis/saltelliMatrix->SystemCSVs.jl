import Pkg
Pkg.activate("../")


using Revise
import PartisanLearning as pl
 const pla = pl.PartyLabel
import CSV

designmatrix = CSV.read("../../../data/saltelli_matrix.csv",
                        pla.DF.DataFrame)

simulation_constants = Dict(:nagents => 500, :nissues => 2, :niterations => 1000)


pla.run_analysis(simulation_constants, designmatrix, 1)
