import Pkg
Pkg.activate("../")


using Revise
import PartisanLearning as pl
const pla = pl.PartyLabel
JULIA_PYTHONCALL_EXE = "/home/marcelovmaciel/miniconda3/bin/python"

using PythonCall

# PythonCall.Deps.add(conda_channels = ["conda-forge"],
#                     conda_packages = ["SALib"]
#                     )

varnames =  ["ncandidates",  "κ", "δ"]
bounds = [[2.,20.], [0.,20.], [1.,10.] ]


function boundsdict_toparamsdf(varnames, bounds;samplesize = 5_000)
    saltelli = pyimport("SALib.sample.saltelli")
    function saltellidict(varnames::Vector{String}, bounds)
        Dict("names" => varnames,
             "bounds" => PythonCall.PyList(bounds),
             "num_vars" => length(varnames))
    end

    boundsdict = saltellidict(varnames, bounds)

    problem_array = pyconvert(Array,
                              saltelli.sample(boundsdict, samplesize,
                                              calc_second_order = true ))


    foodf = pla.DF.DataFrame()

    for (index,value) in enumerate(boundsdict["names"])
        setproperty!(foodf,Symbol(value), problem_array[:, index])
    end
    return(foodf)
end

boundsdict_toparamsdf(varnames, bounds)


n = 5000
lb = [first(param) for param in [ncandidates, κ, δ]]
ub = [param[2] for param in [ncandidates, κ, δ]]


As = QuasiMonteCarlo.sample(n,lb,ub,SobolSample())

Int.(round.(As[1, :]))

Pkg.build("PyCall")
