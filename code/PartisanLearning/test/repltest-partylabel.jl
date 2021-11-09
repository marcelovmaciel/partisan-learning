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

# * PartyLabelModel
# **  BUG: run the δ thing might be leading to a bug!
ncandidates = 10
nissues = 2

m = pla.initialize_model(1000,nissues,
                         ncandidates, δ = 2)

pla.DummyStreakCounter()

([Distances.euclidean(m[k].pos,
                     m[m.properties[:partiesposs][k][:partycandidate]].pos)
     for k in keys(m.properties[:partiesposs])])

pla.model_step!(m)

maximum([Distances.euclidean(m[k].pos,
                     m[m.properties[:partiesposs][k][:partycandidate]].pos)
     for k in keys(m.properties[:partiesposs])])

# ** Try to analyze
ncandidates = 10
nissues = 2

m.properties


m = pla.initialize_model(1000,nissues,
                         ncandidates, δ = 3.)
params = Dict(:κ => 0.0:1.:20.,
              :δ => 3.:1.:20.)

agent_colors(a) = a.id == m.properties[:incumbent]  ? :yellow : (a.amIaCandidate  ?  "#bf2642"  : "#2b2b33")
agent_size(a) = a.id == m.properties[:incumbent]  ? 20 : (a.amIaCandidate ? 15 : 5)


#higher kappa means agents will tend to vote more for
#their partyid. Conversely, lower kappa means people vote more
# for candidates closer to then rather than closer to their party
# higher δ means parties sample further from their location.
# both depended upon the underlying boundaries! think about that !!!

adata = [(a->(pla.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
         (a->(pla.get_distance_IvsPartyCandidate(a,m)), d -> pla.get_representativeness(d,m))]

alabels = ["Against PartyId", "Representativeness"]
mlabels = ["Incumbent", "ENP", "Streaks"]

fig,adf,mdf = abm_data_exploration(m,
                                   pla.abm.dummystep,
                                   pla.model_step!,
                                   params;
                                   adata, pla.mdata,
                                   alabels,mlabels,
                                   ac = agent_colors,
                                   as = agent_size, spu = 1
                                   , static_preplot! = pla.static_preplot! )


pla.get_ENP(m)


# ** Other tests

ncandidates = 10
nissues = 5


m = pla.initialize_model(1000,nissues,
                         ncandidates, δ = 1)

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
