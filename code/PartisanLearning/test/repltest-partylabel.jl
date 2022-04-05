import Pkg

Pkg.activate("../")

using Revise
import PartisanLearning as pl
const is = pl.IssueSalience
const pid = pl.PartyId
const pl.= pl.PartyLabel
using GLMakie
using Agents
import Distances
using CSV
using DataFrames
import ColorSchemes
using StatsBase
include("../test/visualize_model.jl")


# ** Try to analyze

ncandidates = 2
nissues = 2


# FIXME get_plurality_result is giving
# me errors of empty collection in dictmap
# error in setting the candidates maybe
# gotta debug this tomorrow. I'm so tired of this shit




m = pl.initialize_model(1000,
                         nissues,
                         ncandidates,
                         δ=15,
                         κ = 0.,
                         switch=:plurality,
                         ω = 0.99,
                         kappa_switch= :off,
                         special_bounds = (true, (100., 5.)),
                         voter_pos_initializor = () ->
                             pl.sample_1dnormal((100., 5.),
                                                 pl.standard_1d_poss ),
                         party_pos_hardwired = false)


m.properties[:party_switches]


#visualize_noslider(m)

visualize_simpler(m)

for _ in 1:100
    pl.model_step!(m)
end


pl.dist.euclidean((0,2), (100,2))


pl.dictmap(pl.proportionmap,(m.properties[:voterBallotTracker]))
m.properties[:withinpartyshares]



# FIXME: doesnt make sense given the voterballot tracker
# so maybe im not updating the withinpartyshares correctly and this leads
# to some weird behavior downstream

# Save data for quick plotting

for i in pl.abm.allids(m)
    m[i].myPartyId |> println
end


m.properties[:voterBallotTracker]

m.properties[:parties_ids]

m.properties[:voters_partyids]

m.properties




visualize_model(m)

m.properties[:κ]

pl.model_step!(m)


m[2]

pl.dictmap(v->pl.get_mean_among_supporters(v,m),
            pl.get_parties_supporters(m))

# ** Other tests

ncandidates = 10
nissues = 5


m = pl.initialize_model(1000,nissues,
                         ncandidates, δ = 1)




m[1] |> typeof |> fi

pl.model_step!(m)


pl.candidates_iteration_setup!(m)

proportion_peers_voteLikeMe2(1,m)

withinpartyshares[m.properties[:voters_partyids][1]]



m.properties[:voterBallotTracker]

filter((t->t[2]==810),collect(m.properties[:voters_partyids]))

collect(m.properties[:voters_partyids])

m.properties[:voters_partyids]

# *** Test plotting
ncandidates = 5
nissues = 2


m = pl.initialize_model(1000,nissues,
                         ncandidates, δ = 1, κ=0.8)


agent_colors(a) = a.id == m.properties[:incumbent]  ? :yellow : (a.amIaCandidate  ?  "#bf2642"  : "#2b2b33")
agent_size(a) = a.id == m.properties[:incumbent]  ? 20 : (a.amIaCandidate ? 15 : 5)

adata = [(a->(pl.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
         (a->(pl.get_distance_IvsParty(a,m)), pl.StatsBase.mean)]

alabels = ["Voted Against PartyId", "dist(closest,party's)"]

fig,adf,mdf = abm_data_exploration(m,
                                   pl.abm.dummystep,
                                   pl.model_step!,
                                   Dict();
                                   adata, pl.mdata,
                                   alabels,
                                   ac = agent_colors,
                                   as = agent_size, spu = 1)




# * Hardwired case

ncandidates = 2
nissues = 2



m = pl.initialize_model(1000,nissues, ncandidates, δ=5, κ = 28. , switch
    =:runoff, ω = 0.8, kappa_switch= :off,special_bounds = (true, (100., 5.)),
                         voter_pos_initializor = () ->
                             pl.sample_1dnormal((100., 5.),
                                                 pl.one_modal_dispersed ),
                         party_pos_hardwired = true)

loyalty(i) = pl.get_keep_party_id_prob(i.id,m)
ideal_point(i) = i.pos[1]

adata = [ideal_point, loyalty, :myPartyId]

for i in 1:4
    pl.model_step!(m)
end

df,data_m= run!(m, pl.abm.dummystep, pl.model_step!, 1000; adata)
# CSV.write("../../../data/hardwired.csv", data_a)
#df = CSV.read("../../../data/hardwired.csv", DataFrame)

f
f = Figure()
ax = Axis(f[1, 1], xlabel = "step", ylabel = "Party Loyalty")


for i in 1:100
    subsetted_df = filter(:id => x->x==i,df )
    color  = get(ColorSchemes.thermal,
                 (subsetted_df.ideal_point |> unique |> first)/100 )
    scatter!(subsetted_df.step,
             subsetted_df.loyalty, color = (color, 0.2),
               axis = (aspect = 1, xlabel = "x axis", ylabel = "y axis"))

end


Colorbar(f[1, 2], limits = (0, 1), colormap = ColorSchemes.thermal,
         label = "Ideal Point Value")
f


save("first_run_k24_rho15.png", f)





# * Visualize things he asks for

ncandidates = 2
nissues = 2



m = pl.initialize_model(1000,
                         nissues,
                         ncandidates,
                         δ=30,
                         κ = 27.,
                         switch=:plurality,
                         ω = 0.8,
                         kappa_switch= :off,
                         special_bounds = (true, (100., 5.)),
                         voter_pos_initializor = () ->
                             pl.sample_1dnormal((100., 5.),
                                                 pl.more_dispersed_1d_poss),
                         party_pos_hardwired = false)

loyalty(i) = pl.get_keep_party_id_prob(i.id,m)
ideal_point(i) = i.pos[1]

adata = [ideal_point, loyalty, :myPartyId]

for i in 1:4
    pl.model_step!(m)
end

df,data_m= run!(m, pl.abm.dummystep, pl.model_step!, 20; adata)
# CSV.write("../../../data/hardwired.csv", data_a)
#df = CSV.read("../../../data/hardwired.csv", DataFrame)


f = Figure()
ax = Axis(f[1, 1], xlabel = "step", ylabel = "Party Loyalty")


for i in 1:100
    subsetted_df = filter(:id => x->x==i,df )
    color  = get(ColorSchemes.thermal,
                 (subsetted_df.ideal_point |> unique |> first)/100 )
    scatter!(subsetted_df.step,
             log.(subsetted_df.loyalty ), color = (color, 0.4),
               axis = (aspect = 1, xlabel = "x axis", ylabel = "y axis"))

end



Colorbar(f[1, 2], limits = (0, 1), colormap = ColorSchemes.thermal,
         label = "Ideal Point Value")
f


save("log_hardwired_turnout_kappa27.png", f)


