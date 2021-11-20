using Distributed
addprocs(2)


@everywhere import Pkg
@everywhere Pkg.activate("../")
@everywhere JULIA_PYTHONCALL_EXE = "/home/marcelovmaciel/miniconda3/bin/python"

@everywhere using Revise
@everywhere import PartisanLearning as pl
@everywhere const pla = pl.PartyLabel



# PythonCall.Deps.add(conda_channels = ["conda-forge"],
#                     conda_packages = ["SALib"]
#                     )

@everywhere varnames =  ["ncandidates",  "κ", "δ"]
@everywhere bounds = [[2.,20.], [0.,20.], [1.,10.] ]

@everywhere design_matrix = pla.boundsdict_toparamsdf(varnames, bounds)

@everywhere design_matrix[!, :ncandidates] = Int.(round.(design_matrix[!, :ncandidates]))

design_matrix

@everywhere nagents = 1000
@everywhere nissues = 2
@everywhere niterations = 3000

@everywhere adata = [(a->(pla.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
         (a->(pla.get_distance_IvsPartyCandidate(a,m)), d -> pla.get_representativeness(d,m))]

@everywhere mdata = [pla.normalized_ENP,
         x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
         x-> x.properties[:party_switches][end]/x.properties[:nagents]]

for i in 1:10
@everywhere x = eachrow(design_matrix)[1]
    ms = [pla.initialize_model( nagents, nissues, x.ncandidates,
                               κ=  x.κ, δ =  x.δ,
                               seed = s) for s in rand(UInt8, 3)]
    pla.abm.ensemblerun!(ms, pla.abm.dummystep, pla.model_step!, 5;
                     adata,
                     mdata, parallel = true )
end
