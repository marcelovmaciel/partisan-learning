import Pkg

Pkg.activate("../")

using Revise
import PartisanLearning as pl
const is = pl.IssueSalience
const pid = pl.PartyId
const pla = pl.PartyLabel
using GLMakie
using Agents
import Distances
using CSV
using DataFrames
import ColorSchemes
using StatsBase
include("../test/visualize_model.jl")


function get_data_initial_dist(whichdist, nsteps= 20)
    ncandidates = 2
    nissues = 2
    m = pla.initialize_model(1000,
                             nissues,
                             ncandidates,
                             δ=30,
                             κ = 0.,
                             switch=:runoff,
                             ω = 0.8,
                             kappa_switch= :off,
                             special_bounds = (true, (100., 5.)),
                             voter_pos_initializor = () ->
                                 pla.sample_1dnormal((100., 5.),
                                                     whichdist),
                             party_pos_hardwired = false)

    # #(a->(pla.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]
    #

    for _ in 1:4
        pla.model_step!(m)
    end

    loyalty(i) = pla.get_keep_party_id_prob(i.id,m)
    ideal_point(i) = i.pos[1]

    function get_prop_voted_against(m)
        [pla.HaveIVotedAgainstMyParty(a,m)
         for a in allids(m)] |>
             (x-> count(x)/m.properties[:nagents])
    end


    adata = [ideal_point,
             loyalty]

    mdata = [get_prop_voted_against]

    _,data_m= run!(m,
                    pla.abm.dummystep,
                    pla.model_step!, nsteps;
                    adata,
                    mdata)
    return(data_m)
end


function plot_scatter(df, scattername, iter)
    f2 = Figure()
    ax = Axis(f2[1, 1],
          xlabel = "step",
              ylabel = "Prop Cross Party Vote")
    scatter!(df.step,
             df.get_prop_voted_against)
    stringtosave = "./plots/" * scattername * string(iter) * ".png"
    save(stringtosave, f2)
end


function summaries_prop_voted_against(df)
    DataFrame(Dict(zip(
        [:min,
         :max,
         :mean],
        [minimum(df[:get_prop_voted_against]),
         maximum(df[:get_prop_voted_against]),
         StatsBase.mean(df[:get_prop_voted_against])])))
end


function collect_prop_voted_against(whichdist)

    holderdf = DataFrame([Int64[],
                          Float64[],
                          Int64[]],
                         [:step,
                          :get_prop_voted_against,
                          :iter])
    for i in 1:30
        df = get_data_initial_dist(whichdist)
        df[:iter] = [i for _ in 1:21]

        holderdf = vcat(holderdf,
                        df)
        #plot_scatter(df, scattername, i)
    end

    return(holderdf)
end

overlap20df = collect_prop_voted_against(pla.overlap_20_poss)
overlap50df = collect_prop_voted_against(pla.standard_1d_poss)
overlap80df = collect_prop_voted_against(pla.overlap_80_poss)



f2 = Figure()
ax = Axis(f2[1, 1], xlabel = "step",
          ylabel = "Prop Cross Party Vote")


boxplot!(overlap80df.step,
             overlap80df.get_prop_voted_against)
