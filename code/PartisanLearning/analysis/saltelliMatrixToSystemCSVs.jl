import Pkg
Pkg.activate("../")

using Revise
import PartisanLearning as pl
 const pla = pl.PartyLabel
import CSV

designmatrix = CSV.read("../../../data/saltelli_matrix.csv",
                        pla.DF.DataFrame)
simulation_constants = Dict(:nagents => 5000,
                            :nissues => 2,
                            :niterations => 200)

#test
pla.run_analysis(simulation_constants, designmatrix, 1)
