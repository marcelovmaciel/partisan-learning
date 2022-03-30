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

m_params = pla.ModelParams()
m = pla.initialize_model(m_params)

# * Vis functions
foursteps!(m) =  for _ in 1:4 pla.model_step!(m) end

function mean_loyalty(m)
[pla.get_keep_party_id_prob(i,m) for i in pla.abm.allids(m)] |>  pla.mean
        end

function prop_pswitch(m)
    m.properties[:party_switches][end]/m.properties[:nagents]
 end

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
    #           :kappa_switch => [:off, :on])

function single_interactive_vis(m)

    foursteps!(m)
    adata = [(a->(pla.get_distance_IvsPartyCandidate(a,m)), d -> pla.get_representativeness(d,m))]

    mdata =  [pla.get_2candidates_distance, prop_pswitch, mean_loyalty]

    alabels = ["rep"]
    mlabels = ["dist c1c2", "prop_switch","mean_loyalty"]
    fig,pl2 = abmexploration(m;
                             agent_step! = pla.abm.dummystep,
                             model_step! = pla.model_step!,
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

fig
