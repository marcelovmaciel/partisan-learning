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
using AlgebraOfGraphics

include("../test/visualize_model.jl")


function get_data_initial_dist(whichdist, nsteps= 20)
    ncandidates = 2
    nissues = 2
    m = pla.initialize_model(1000,
                             nissues,
                             ncandidates,
                             δ=15,
                             κ = 0.1,
                             switch=:plurality,
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

    function prop_pswitch(m)

             m.properties[:party_switches][end]/m.properties[:nagents]
    end

    function prop_crossvoting(m)
      m.properties[:cross_voting][end]/m.properties[:nagents]
    end

    function mean_loyalty(m)
[pla.get_keep_party_id_prob(i,m) for i in pla.abm.allids(m)] |>  pla.mean
        end


    adata = [ideal_point,
             loyalty]

    mdata = [prop_pswitch, mean_loyalty, prop_crossvoting]

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


# function summaries_prop_voted_against(df)
#     DataFrame(Dict(zip(
#         [:min,
#          :max,
#          :mean],
#         [minimum(df[:get_prop_voted_against]),
#          maximum(df[:get_prop_voted_against]),
#          StatsBase.mean(df[:get_prop_voted_against])])))
# end


function collect_prop_voted_against(whichdist)

    holderdf = DataFrame([Int64[],
                          Float64[],
                          Int64[], Float64[], Float64[]],
                         [:step,
                          :prop_pswitch,
                          :iter, :mean_loyalty, :prop_crossvoting])
    for i in 1:100
        df = get_data_initial_dist(whichdist)
        df[!,:iter] = [i for _ in 1:21]

        holderdf = vcat(holderdf,
                        df)
        #plot_scatter(df, scattername, i)
        #println("\n", "new iteration\n")
    end

    return(holderdf)
end


overlap20df = collect_prop_voted_against(pla.overlap_20_poss)

overlap50df = collect_prop_voted_against(pla.standard_1d_poss)

overlap80df = collect_prop_voted_against(pla.overlap_80_poss)


overlap20df[!,:type] =[:o20 for _ in 1:(overlap20df |> size)[1]]
overlap50df[!,:type] =[:o50 for _ in 1:(overlap50df |> size)[1]]
overlap80df[!,:type] = [:o80 for _ in 1:(overlap80df |> size)[1]]


gdf = vcat(overlap20df, overlap50df, overlap80df)



plt = data(gdf) * mapping(:step,
                          :prop_pswitch,
                          color = :type => nonnumeric,
                          dodge = :type => nonnumeric) * visual(BoxPlot);

foo = draw(plt)


stringtosave = "./plots/overlapping_pswitch.png"

AlgebraOfGraphics.save(stringtosave, foo, px_per_unit = 5)


plt2 = data(gdf) * mapping(:step,
                          :prop_crossvoting,
                          color = :type => nonnumeric,
                          dodge = :type => nonnumeric) * visual(BoxPlot);

foo2 = draw(plt2)



stringtosave = "./plots/overlapping_test_cross.png"

AlgebraOfGraphics.save(stringtosave, foo2, px_per_unit = 5)

# f2 = Figure()
# ax = Axis(f2[1, 1], xlabel = "step",
#           ylabel = "Prop Cross Party Vote")


# b1 = boxplot!(overlap20df.step,
#              overlap20df.get_prop_voted_against)
# b2 = boxplot!(overlap50df.step,
#               overlap50df.get_prop_voted_against)
# b3 = boxplot!(overlap80df.step,
#               overlap80df.get_prop_voted_against)


# Legend(f2[1,2], [b1,b2,b3], ["20", "50", "80"])




# f3 = Figure()
# ax = Axis(f3[1, 1], xlabel = "step",
#           ylabel = "Prop Cross Party Vote")


# b = boxplot!(gdf.step,
#              gdf.get_prop_voted_against, dodge = gdf.type, color = gdf.type)



# Legend(f3[1,2], [b1,b2,b3], ["20", "50", "80"])
