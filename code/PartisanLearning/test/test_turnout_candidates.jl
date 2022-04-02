import Pkg
Pkg.activate("../")

using Revise
import PartisanLearning as pl


using GLMakie
using Agents
import Distances
using CSV
using DataFrames
import ColorSchemes
using StatsBase
using AlgebraOfGraphics
using ProgressMeter



# include("../test/visualize_model.jl")

m_params = pl.ModelParams()


m = pl.initialize_model(m_params)

# * Vis functions
p1 = m.properties[:parties_ids][1]
p2 = m.properties[:parties_ids][2]

agent_size(a) = (a.id in m.properties[:parties_ids] ?  37 : (a.amIaCandidate ? 35 : 5))
agent_marker(a) = if a.id in m.properties[:parties_ids] '♠' else '∘' end

function agent_colors_2parties(a)
    if a.id == m.properties[:incumbent_party]
        :yellow
    elseif a.myPartyId == p1
        :red
    elseif a.myPartyId == p2
        :blue
    end
end

# params = Dict(:κ => 0.0:1:100.0,
    #           :switch => [:random, :plurality, :runoff],
    #           :δ=> 0.0:5:70,
    #           :ω => 0.0:0.1:1.0,
    #           :kappa_switch => [:off, :on]) FIXME: this breaks the code


function single_interactive_vis(m)
    pl.foursteps!(m)
    adata = [(a->(pl.get_distance_IvsPartyCandidate(a,m)), d -> pl.get_representativeness(d,m))]

    mdata =  [pl.get_2candidates_distance, pl.prop_pswitch, pl.mean_loyalty]

    alabels = ["rep"]
    mlabels = ["dist c1c2", "prop_switch","mean_loyalty"]
    fig,pl2 = pl.abmexploration(m;
                             agent_step! = pl.abm.dummystep,
                             model_step! = pl.model_step!,
                            # adata, # BUG: this breaks the code !
                              mdata,
                             # alabels,
                              mlabels,
                             ac = agent_colors_2parties,
                             as = agent_size,
                             am = agent_marker,
                             spu = 1:100)

    return(fig, pl2)
end


# * Actual test

fig,f = single_interactive_vis(m)



function get_data_initial_dist(params, nsteps= 20)

    m = pl.initialize_model(params)

    pl.foursteps!(m)


    adata = [pl.f_ideal_point]# ,
             # loyalty] FIXME: this breaks the code!

    mdata = [pl.prop_pswitch,
             pl.mean_loyalty,
             pl.prop_crossvoting,
             pl.get_2candidates_distance]

    _,data_m= run!(m,
                    pl.abm.dummystep,
                    pl.model_step!, nsteps;
                    adata,
                    mdata)
    return(data_m)
end



function collect_runs(params)

    holderdf = DataFrame([Int64[],
                          Float64[],
                          Int64[],
                          Float64[],
                          Float64[],
                          Float64[]],
                         [:step,
                          :prop_pswitch,
                          :iter,
                          :mean_loyalty,
                          :prop_crossvoting,
                          :get_2candidates_distance])

@showprogress 1 "Running "  for i in 1:100
        df = get_data_initial_dist(params)
        df[!,:iter] = [i for _ in 1:21]

        holderdf = vcat(holderdf,
                        df)
    end

    return(holderdf)
end


function collect_per_overlap(whichdist)
    (whichdist |>
        pl.overlap_initializor |>
        x-> pl.ModelParams(voter_pos_initializor = x) |>
        collect_runs)
end

overlap20df = collect_per_overlap(pl.overlap_20_poss)

overlap50df = collect_per_overlap(pl.overlap_50_poss)

overlap80df = collect_per_overlap(pl.overlap_80_poss)


function concat_overlap_data(f,s,t)

    for (i,j) in zip([f,s,t], [:o20, :o50, :o80])
        println(j)
        @show first(i)

        df_size =  (i |> size)[1]
        i[!,:type] =[j for _ in 1:df_size]
    end
    gdf = vcat(f,s,t)
    return(gdf)
end


#m_params = pl.ModelParams(\k)
#m_params |> typeof |> fieldnames



gdf =  concat_overlap_data(overlap20df,
                           overlap50df,
                           overlap80df)

plt = data(gdf) * mapping(:step,
                          :get_2candidates_distance,
                          color = :type => nonnumeric,
                          dodge = :type => nonnumeric) * visual(BoxPlot);

foo = draw(plt)

stringtosave = "./plots/overlapping_2cdist_kappa1.png"

AlgebraOfGraphics.save(stringtosave, foo, px_per_unit = 5)
