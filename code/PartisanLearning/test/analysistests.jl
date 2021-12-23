import Pkg
Pkg.activate("../")
JULIA_PYTHONCALL_EXE = "/home/marcelovmaciel/miniconda3/bin/python"



using Revise
import PartisanLearning as pl
const pla = pl.PartyLabel
import CSV
using Base.Filesystem
using GLMakie
using PythonCall

sobol = pyimport("SALib.analyze.sobol")
PythonCall.Deps.add(conda_channels = ["conda-forge"],
                      conda_packages = ["jupyterlab"]
                    )

varnames =  ["ncandidates",  "κ", "δ"]
bounds = [[2.,15.], [0.,7.], [0.5,7.]]
testdict = pla.saltellidict(varnames,bounds)

outputnames = ["Eccentricity",
                   "NENP",
                   "PartySwitches",
                   "Representativeness", "LongestIStreak"]

designmatrix = CSV.read("../../../data/saltelli_matrix.csv",
                        pla.DF.DataFrame)


outputmatrix = CSV.read("../../../data/output_matrix.csv",
                        pla.DF.DataFrame)

PythonCall.Py(outputmatrix[!,outputnames[1]])

sobol.analyze(testdict,
              outputmatrix[!,outputnames[1]])

PythonCall.PyArray(outputmatrix[!,outputnames[1]], array = true).size

sio = hcat(designmatrix,outputmatrix)


designmatrix
outputmatrix

function scatter_grid(inputvar,whichf, data = sio)
    f = Figure()
    outputnames = ["Eccentricity",
                   "Representativeness",
                   "NENP",
                   "PartySwitches",
                    "LongestIStreak"]
    poss = [(1,1), (1,2), (1,3), (2,1), (2,2:3)]
    for (outputname,pos) in zip(outputnames, poss)
        GLMakie.scatter!(Axis(whichf[pos...], xlabel = inputvar, ylabel = outputname ),
                         data[!,inputvar],data[!, outputname])
    end

end

f = Figure()
scatter_grid("δ", f)

GLMakie.save("δ_out.png", f)

f = Figure()
scatter_grid("κ", f)
GLMakie.save("κ_out.png", f)

f = Figure()
scatter_grid("ncandidates", f )
GLMakie.save("ncandidates_out.png", f)
