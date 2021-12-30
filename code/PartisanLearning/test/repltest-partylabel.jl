import Pkg

Pkg.activate("../")


using Revise
import PartisanLearning as pl
const is = pl.IssueSalience
const pid = pl.PartyId
const pla = pl.PartyLabel
using InteractiveDynamics
using GLMakie
using Agents
import Distances



# ** Try to analyze
ncandidates = 5
nissues = 2

m = pla.initialize_model(500,nissues,
                         ncandidates, δ = 3.)

foo = pla.select_primariesCandidates(m.properties[:partiesposs],m)

pla.get_primaries_votes(m,foo)

foo2 = pla.dictmap(pla.proportionmap,pla.get_primaries_votes(m,foo))



get_plurality_result(pla.get_primaries_votes(m,foo))
get_runoff_result(pla.get_primaries_votes(m,foo),m)



#Remember: INCUMBENT IS THE CANDIDATE!!!!!!!

delete!(parties_supporters, m[m.properties[:incumbent]].myPartyId)
pla.dictmap(parties_supporters,
        supporters->map( i ->pla.get_closest_fromList(i,
                                                      foo[m[i].myPartyId], m)
                         ,supporters))

map(i->pla.get_closest_fromList(i, foo[m[i].myPartyId], m), parties_supporters[243])

parties_supporters

foo

m.properties[:partiesposs]

pla.dictmap(parties_supporters,
        supporters->map( i ->pla.get_closest_fromList(i,
                                                      foo[m[i].myPartyId], m)
                         ,supporters))


parties[1]

#=
- I sweep over agents
- Agents tell which one is closer to them
- I keep track of that
- I calculate the proportions
=#

m.properties

m.properties[:voters_partyids]

foo |> typeof

params = Dict(:κ => 0.0:1.:7,
              :δ => 0.5:0.5:7)


agent_colors(a) = a.id == m.properties[:incumbent]  ? :yellow : (a.amIaCandidate  ?  "#bf2642"  : "#2b2b33")
agent_size(a) = a.id == m.properties[:incumbent]  ? 20 : (a.amIaCandidate ? 15 : 5)


#higher kappa means agents will tend to vote more for
#their partyid. Conversely, lower kappa means people vote more
# for candidates closer to then rather than closer to their party
# higher δ means parties sample further from their location.
# both depended upon the underlying boundaries! think about that !!!

adata = [#(a->(pla.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
         (a->(pla.get_distance_IvsPartyCandidate(a,m)), d -> pla.get_representativeness(d,m))]

mdata = [pla.normalized_ENP,
         x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
         x-> x.properties[:party_switches][end]/x.properties[:nagents],
         pla.get_incumbent_eccentricity]

alabels = ["Rep"]
mlabels = ["NENP", "IStreaks", "PSwitches", "ecc"]

fig,adf,mdf = abm_data_exploration(m,
                                   pla.abm.dummystep,
                                   pla.model_step!,
                                   params;
                                   adata, mdata,
                                   alabels,mlabels,
                                   ac = agent_colors,
                                   as = agent_size, spu = 1 )

# TODO: think about this weird behavior in which people switch
# more than they vote against their party
# Maybe I'm doing some mistake on counting the switches



function pltfn(foo, fig = fig)

k =  collect(keys(foo.val))
pie(fig[1,3], [foo.val[v] for v in k], color = k, radius = 5, inner_radius = 2)
             end




partyproportions = pla.proportionmap([v for (k,v) in m.properties[:voters_partyids]])
partyids = collect(keys(partyproportions))
foo = Observable([partyproportions[p] for p in partyids])


function newstep(m, foo = foo)
    pla.model_step!(m)
    partyproportions = pla.proportionmap([v for (k,v) in m.properties[:voters_partyids]])
    foo[] = [partyproportions[p] for p in partyids]
end



fig,adf,mdf = abm_data_exploration(m,
                                   pla.abm.dummystep,
                                   newstep,
                                   params;
                                   adata, mdata,
                                   alabels,mlabels,
                                   ac = agent_colors,
                                   as = agent_size, spu = 1 )

barplot(fig[1,3], 1:length(partyids), foo)


# ** Other tests

ncandidates = 10
nissues = 5


m = pla.initialize_model(1000,nissues,
                         ncandidates, δ = 1)

foreach(Pkg.add, (  "GLMakie", "InteractiveDynamics"))



m[1] |> typeof |> fi

pla.model_step!(m)


pla.candidates_iteration_setup!(m)

proportion_peers_voteLikeMe2(1,m)

withinpartyshares[m.properties[:voters_partyids][1]]



m.properties[:voterBallotTracker]

filter((t->t[2]==810),collect(m.properties[:voters_partyids]))

collect(m.properties[:voters_partyids])

m.properties[:voters_partyids]

# *** Test plotting
ncandidates = 5
nissues = 2


m = pla.initialize_model(1000,nissues,
                         ncandidates, δ = 1, κ=0.8)


agent_colors(a) = a.id == m.properties[:incumbent]  ? :yellow : (a.amIaCandidate  ?  "#bf2642"  : "#2b2b33")
agent_size(a) = a.id == m.properties[:incumbent]  ? 20 : (a.amIaCandidate ? 15 : 5)

adata = [(a->(pla.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
         (a->(pla.get_distance_IvsParty(a,m)), pla.StatsBase.mean)]

alabels = ["Voted Against PartyId", "dist(closest,party's)"]

fig,adf,mdf = abm_data_exploration(m,
                                   pla.abm.dummystep,
                                   pla.model_step!,
                                   Dict();
                                   adata, pla.mdata,
                                   alabels,
                                   ac = agent_colors,
                                   as = agent_size, spu = 1)


# ** Eccentricity test

ncandidates = 15
nissues = 2

m = pla.initialize_model(500,nissues,
                         ncandidates, δ = 3.)
pla.get_median_pos(m)

dims = m.properties[:nissues]
medians = Vector{Float64}()

for dim in 1:dims
    dimmedian = Statistics.median([m[x].pos[dim] for x in  pla.abm.allids(m)])
        push!(medians, dimmedian)
end


# TODO: check Incumbent again


m.properties[:median_pos]
pla.dist.euclidean(m.properties[:incumbent])


# * Test sampler

ncandidates = 2
nissues = 2

m = pla.initialize_model(500,nissues,
                         ncandidates, δ = 3.)
params = Dict(:κ => 0.0:1.:20.,
              :δ => 1:1:20)

# TODO: run with 0.5,0.5 below : wtf is that a periodic boundary condition on?
# FIXME: yes, periodic conditions are on! LOL
neighs = pla.abm.nearby_ids((10,10),
                       m,
                       1,
                       exact = true) |> collect

[pla.dist.euclidean((5,5), m[x].pos) for x in neighs] |> maximum

agent_colors(a) = a.id in neighs  ? :yellow : "#2b2b33"
agent_size(a) = a.id in neighs  ? 20 :  5


#higher kappa means agents will tend to vote more for
#their partyid. Conversely, lower kappa means people vote more
# for candidates closer to then rather than closer to their party
# higher δ means parties sample further from their location.
# both depended upon the underlying boundaries! think about that !!!

adata = [(a->(pla.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
         (a->(pla.get_distance_IvsPartyCandidate(a,m)), d -> pla.get_representativeness(d,m))]

mdata = [pla.normalized_ENP,
         x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
         x-> x.properties[:party_switches][end]/x.properties[:nagents],
         pla.get_incumbent_eccentricity]

alabels = ["¬-PartyId", "Rep"]
mlabels = ["NENP", "IStreaks", "PSwitches", "ecc"]

fig,adf,mdf = abm_data_exploration(m,
                                   pla.abm.dummystep,
                                   pla.model_step!,
                                   params;
                                   adata, mdata,
                                   alabels,mlabels,
                                   ac = agent_colors,
                                   as = agent_size, spu = 1
                                   , static_preplot! = pla.static_preplot! )
