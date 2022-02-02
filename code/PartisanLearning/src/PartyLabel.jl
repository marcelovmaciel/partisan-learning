# * Party Label model

#= This is a copy of the Party id model. Why am I copying
it? Because the previous version is not wrong, but design needs to change. I
still wanna play with the other version so I'll keep it. =#
#=
What will change in this version?
I'll need to add some sticky parties,
rather than having parties as a latent thing
in the model. So this changes the initial condition.
I'll also need to change how people set their party position...
Also, the update rule is completely different now.
See the notes/sketchsofaModel/partyid-sketch.pdf design
document for more on that.
=#
module PartyLabel
# ** Initial condition
#= the initial logic is the following:
• One initializes the voters;
• Some c number of voters will be treated as “candidates”;
• Each voter votes for the candidate who is closest to them. The one
with the most votes becomes the incumbent;
• Maybe each voter treats their candidate id as their initial partyid?
This is a model initialization artifact=#

import Agents as abm
import Distributions as distri
import Distances as dist
import Base.@kwdef
using StaticArrays
using StatsBase
import Statistics
# import GLMakie
import DataFrames as DF
using PythonCall
using Random
import CSV
using ProgressMeter
using Distributions

mutable struct Voter{n} <: abm.AbstractAgent
    id::Int
    pos::NTuple{n,Float64}
    amIaCandidate::Bool
    κ::Float64
    myPartyId::Int
end

  #= Note that if a voter is a candidate then its ~myPartyId~ should be the
agent's id. Maybe I'll create an dictionary Int=> Symbol to identify the parties
throughout simulation inspection =#
const bounds = (0,100)

@kwdef struct ModelParams
    nagents = 3000
    ncandidates = 2
    nissues = 2
    bounds = bounds
    κ = 1. # this influences whether I'll vote against my party or not
    δ = 3. # this influences the radius of candidates a party samples from
end


function sample_uniform_pos(nissues = 2)
Tuple(rand(distri.Uniform(bounds...),nissues))
end


function sample_1dnormal(bound)
    #=
    - I'll simply put a constant in one dimension and work on the other lol.
- I'll use cohen d to define overlapping distributions:
  - As shown here [[https://rpsychologist.com/cohend/][Interpreting Cohen&#x27;s d | R Psychologist]] a cohen d of 1.35 gives an overlap of 50%.
- Thus if I have a distribution of (43,10) I gotta have the other as (56.5,10)
=#
    if rand([false,true])
        pos = (rand(Normal(43.25,10)), bound[2]/2 )
    else
        pos = (rand(Normal(56.75,10)),bound[2]/2 )
    end
    return(pos)
end


function Voter(id::Int;nissues =  1,
               pos = () -> sample_uniform_pos(nssissues), κ = 10.)
    amIaCandidate = false
    myPartyId = -3
    return(Voter{nissues}(id,pos(),amIaCandidate, κ, myPartyId))
end

 function dictmap(l,d)
    Dict(Pair(k,l(v)) for (k,v) in d)
end

function kvdictmap(l,d)
    Dict(Pair(k,l(k,v)) for (k,v) in d)
end

function kdictmap(l,d)
    Dict(Pair(k,l(k)) for (k,_) in d)
end


function sample_parties_pos(nparties, model)
    ids = collect(abm.allids(model))

    partiesposs = Dict(map(x-> (x,model[x].pos),
                           sample(ids,nparties, replace = false)))
    actualpartiesposs = dictmap(v-> Dict(:partyposition => v,
                                         :partycandidate => -1,
                                         :δ => model.properties[:δ]),
                                partiesposs)
    return(actualpartiesposs)
end


function sample_candidates(party,m,n=4)
    δ = m.properties[:δ]

    nearby_agents = collect(abm.nearby_ids(m[party], m, δ, exact = true))

    nearby_agents_from_the_party = filter(agentid -> m[agentid].myPartyId == party,
                                          nearby_agents)
    (isempty(nearby_agents_from_the_party) ?
        fill(party, (n,1)) :
        sample(nearby_agents_from_the_party,n))
end


function select_primariesCandidates(model::abm.ABM)
    party_candidate_pairs = (
        model.properties[:parties_ids] .|>
            (pid -> Pair(pid,sample_candidates(pid,model))) |>
        Dict)
      return(party_candidate_pairs)
end


function get_plurality_result(primariesresult::Dict)
    dictmap(argmax ∘ proportionmap, primariesresult)
end


function get_plurality_result(m::abm.ABM)
    candidates = select_primariesCandidates(m)
    primariesresult = get_primaries_votes(m,candidates)
    get_plurality_result(primariesresult)
end

function get_runoff_result(primariesresult,m)
    primariesproportion = dictmap(proportionmap
                                  ,primariesresult)
    function runoff(k,v)
             if any(x->x>0.5,values(v))
                 argmax(v)
             else
                 toptwo = sort(collect(v),
                               by=x->x[2],
                               rev = true)[1:2] .|> x->x[1]
                 (map(x->get_closest_fromList(x,toptwo,m),
                  get_parties_supporters(m)[k]) |>
                      proportionmap |>
                      argmax)
             end
    end
    kvdictmap(runoff, primariesproportion)
end


function get_runoff_result(m)
    candidates = select_primariesCandidates(m)
    primariesresult = get_primaries_votes(m,candidates)

    get_runoff_result(primariesresult,m)
end

# FIXME: Should create a radius  + m.properties[:parties_ids] version
function get_random_candidates(m)
    parties_supporters = get_parties_supporters(m)
        if m.properties[:incumbent_party] == 0
            dictmap(rand,parties_supporters)
    else
        delete!(parties_supporters, m[m.properties[:incumbent_party]].myPartyId)
        dictmap(rand,
                parties_supporters)
    end
end

function set_candidates!(model, switch)

    if switch == :runoff
        candidateids = get_runoff_result(model)
    elseif switch == :plurality
        candidateids = get_plurality_result(model)
    else
        candidateids = get_random_candidates(model)
    end


    for (pid, candidateid) in candidateids

        model[candidateid].amIaCandidate = true
        model[candidateid].myPartyId = pid
        model.properties[:parties_candidateid_ppos_δ][pid][:partycandidate] = candidateid
    end

end


"get_closest_candidate(agentid::Int, model::abm.ABM)"
function get_closest_candidate(agentid::Int,model)

    #=
    Those are dummy variables. In the first loop iteration
    they are replaced. Since the bounds: (0,1) there is no
    way distance > dummydistance! the id is negative because
    there can be no such thing. Thus, it will warn me downstream
    if any mistake has happened here.
    =#
    dummydistance = 100000.
    candidateid, itspartyid = -1,-1

    for i in abm.allids(model)     #pid.abm.nearby_ids(model[testid].pos, model)
        if model[i].amIaCandidate
            candidatepos = model[i].pos
             distance = dist.euclidean(candidatepos,
                                      model[agentid].pos)
            if distance < dummydistance
                dummydistance = distance
                candidateid, itspartyid = model[i].id, model[i].myPartyId
            end
        end
    end
    return(candidateid,itspartyid)
end


function get_closest_fromList(agentid,candidate_list,model)
    dummydistance = 100000.
    candidateid = -1
    for candidate in candidate_list
    candidatepos = model[candidate].pos
    distance = dist.euclidean(candidatepos,
                              model[agentid].pos)
        if distance < dummydistance
            dummydistance = distance
            candidateid = candidate
        end
    end
    return(candidateid)
end


function get_primaries_votes(m, primariesCandidatesDict)
    parties_supporters = get_parties_supporters(m)
    get_closest_toI(i) = get_closest_fromList(i,
                         primariesCandidatesDict[m[i].myPartyId],
                                              m)
    return(dictmap(supporters->map(get_closest_toI,supporters),
                parties_supporters))
end

"getmostvoted(model::abm.ABM)"
function getmostvoted(model::abm.ABM, initial_or_iteration = :initial)

    # HACK: this ifelse is for reusing this function in the iteration
    # Maybe later write a tagging type for that
    if initial_or_iteration == :initial
        votingfn = get_closest_candidate
    else
        votingfn = get_whichCandidatePartyAgentVotesfor
    end

    #=The candidate who is closest to most wins.
    Standard downsian assumption.=#
    closest_candidates = Array{Int}(undef, abm.nagents(model))
    for id in abm.allids(model)
        closest_candidates[id] = votingfn(id,model)[1]
    end
   return(argmax(proportionmap(closest_candidates))) # argmax will return only one maximal, beware of that!
end


function get_parties_supporters(model)
    Dict(
        Pair(k,
             collect(keys(filter(t->t[2]==k,
                    model.properties[:voters_partyids]))))
        for k in model.properties[:parties_ids])
end

function get_withinpartyshares(model)
    #= This function is very important and calculates the proportion
    of vote shares WITHIN the parties. That is, to what party
    the supporters of each party have voted for.
    Naturally, it **ONLY MAKES SENSE** if I update
    model.properties[:voterBallotTracker] before using it!!!!=#
    parties_supporters = get_parties_supporters(model)
    withinpartyvotes = dictmap(v-> map(x-> model.properties[:voterBallotTracker][x][end],
                                       v),
                               parties_supporters)
   return(dictmap(proportionmap, withinpartyvotes))
end

mutable struct StreakCounter
    old_incumbentholder::Int64
    has_switchedlist::Vector{Bool}
    current_streak::Int64
    longest_streak::Dict # will hold both streak_value and incumbent_position
end

DummyStreakCounter() = StreakCounter(1,Vector{Bool}(), 0, Dict(:streak_value => 0,
                                                                 :incumbent_pos => (0,)))

function initialize_incumbent_streak_counter!(m)
    m.properties[:incumbent_streak_counter].old_incumbentholder = 0
    m.properties[:incumbent_streak_counter].current_streak = 0
    m.properties[:incumbent_streak_counter].longest_streak[:streak_value] = 0
    m.properties[:incumbent_streak_counter].longest_streak[:incumbent_pos] = ntuple(x->0.,Val(m.properties[:nissues]))
end


function get_median_pos(m)
    dims = m.properties[:nissues]
    medians = Vector{Float64}()
    for dim in 1:dims
        dimmedian = Statistics.median([m[x].pos[dim] for x in  abm.allids(m)])
        push!(medians, dimmedian)
    end
    return(medians)
 end


function initialize_model(nagents::Int, nissues::Int, nparties;
                          κ = 2., δ=1.,  switch= :random,  seed = 125, ω = 0.8, kappa_switch = :off,
                           special_bounds = (false, bounds), voter_pos_initializor = () -> sample_1dnormal(special_bounds[2]))
    if special_bounds[1] == true
        space = abm.ContinuousSpace(special_bounds[2], periodic = false)
    else
        space = abm.ContinuousSpace(ntuple(x -> float(last(bounds)),nissues), periodic = false)
    end


    rng = Random.MersenneTwister(seed)
    # postype = typeof(ntuple(x -> 1.,nissues))

    #=
    There are three main auxiliary model collections:
    - What are the parties ids, and their positions and candidates (partiesposs)
    - What are the voters parties ids NOW (voters_partyids)
    - What are the the voters voting track (votersBallotTracker)
    =#
    voters_partyids = Dict{Int64, Int64}()
    voterBallotTracker = Dict{Int64, Vector{Int}}()
    parties_candidateid_ppos = Dict{}()
    parties_ids = Vector{Int}(undef,nparties)

# This will allow me to stop recalculating it from parties_candidateid_ppos
    withinpartyshares = Dict{}()
    #=
    I am adding that as a model property to later:
    1) use it for synchronous partyid update
    2) plotting (closer to the sugarscape example)
    =#
    properties = Dict(:incumbent_party => 0,
                      :nagents => nagents,
                      :nissues => nissues,
                      :ncandidates => nparties,
                      :parties_ids => parties_ids,
                      :parties_candidateid_ppos_δ => parties_candidateid_ppos,
                      :δ => δ,
                      :switch => switch,
                      :κ => κ,
                      :ω => ω,
                      :voters_partyids => voters_partyids,
                      :voterBallotTracker=>voterBallotTracker,
                      :withinpartyshares => withinpartyshares,
                      :incumbent_streak_counter => DummyStreakCounter(),
                      :keep_probs => [[1.,1]],
                      :party_switches => [0],
                      :median_pos => [],
                      :kappa_switch => kappa_switch,
                      :is_at_step => 0)

    model = abm.ABM(Voter{nissues}, space; rng,properties = properties)

    for i in 1:nagents
        vi = Voter(i, nissues=nissues, κ = model.properties[:κ], pos = voter_pos_initializor)
        abm.add_agent_pos!(vi, model)
    end
    model.properties[:parties_candidateid_ppos_δ] = sample_parties_pos(nparties,
                                                        model)

    for (i,v) in enumerate(collect(keys(model.properties[:parties_candidateid_ppos_δ])))
        model.properties[:parties_ids][i]=v
    end

    for i in abm.allids(model)
        mypartyid = get_closest_fromList(i,model.properties[:parties_ids],model)
        model[i].myPartyId = mypartyid
    end

    model.properties[:voters_partyids] = Dict(
        (model[x].id => model[x].myPartyId)
        for x in abm.allids(model))

    initialize_incumbent_streak_counter!(model)
    model.properties[:voterBallotTracker] = Dict((k,Int64[]) for (k,_) in model.properties[:voters_partyids])
    model.properties[:median_pos] = get_median_pos(model)
    return(model)
end



function assume_initial_partyid!(i, m)


    pairs = collect(proportionmap(m.properties[:voterBallotTracker][719]))
    vals = first.(pairs)
    weights = last.(pairs)
    which_party_id_to_assume = sample(vals, Weights(weights))
    i.myPartyId = which_party_id_to_assume
end

# ** Stepping
#=
there is some sublety here.
#basicamente na iteracao vai ter dois tipos de calculo de distancia pra cada
agente ele quer saber quem ele ta mais proximo em termos de posicao DELE MAS
tambem quem que ta mais proximo da partyid dele entao ele considera dois
candidatos se o candidato mais proximo dele 'e tipo muito mais proximo que o
candidato do partido, uma contante kappa, ele vota no outro
=#
"reset_candidates!(model::abm.ABM)"
function reset_candidates!(model::abm.ABM)
    for agent in abm.allids(model)
        model[agent].amIaCandidate = false
    end
end


"candidates_iteration_setup!(model::abm.ABM)"
function candidates_iteration_setup!(m::abm.ABM)
    reset_candidates!(m)
    set_candidates!(m, m.properties[:switch])
end


function get_whichCandidatePartyAgentVotesfor(agentid, model)
    κ = model.properties[:κ]

    closest_to_me_id_pid = get_closest_candidate(agentid,model)

    mypartycandidate = model.properties[:parties_candidateid_ppos_δ][model[agentid].myPartyId][:partycandidate]

    whoillvotefor = (mypartycandidate, model[agentid].myPartyId)

    two_candidates_distance = dist.euclidean(model[closest_to_me_id_pid[1]].pos,
                                             model[mypartycandidate].pos)


    if two_candidates_distance > κ
        whoillvotefor = closest_to_me_id_pid
    end

    return(whoillvotefor)

end

function get_proportion_peers_voteLikeMe(agentid,model)
    # This function only makes sense after updating m.properties[:voterBallotTracker]
    # and then updating the model.properties[:withinpartyshares]
    shares = model.properties[:withinpartyshares]
    partyImVotingFor = get_whichCandidatePartyAgentVotesfor(agentid,model)[2]
    if !(partyImVotingFor in
         keys(shares[model[agentid].myPartyId]))
        proportion = 0.0
    else
        proportion = shares[model[agentid].myPartyId][partyImVotingFor]
    end
    return(proportion)
end


function HaveIVotedAgainstMyParty(agentid::Int, model)
    mypartycandidate = model.properties[:parties_candidateid_ppos_δ][model[agentid].myPartyId][:partycandidate]
    mycandidate = get_whichCandidatePartyAgentVotesfor(agentid, model)[1]
    return(mycandidate != mypartycandidate)
end

HaveIVotedAgainstMyParty(agent::Voter, model) = HaveIVotedAgainstMyParty(agent.id,model)


#=
1. Pick the proportion of iterations that I voted for the party I’m voting
for in this iteration.
2. Pick the proportion of people from my party that voted
different from me;
- I change my party id to this other party with a probability equal to tanh(proportion 1 + proportion 2). =#

function get_keep_party_id_prob(agentid,model)
    myLast_PartyVote = model.properties[:voterBallotTracker][agentid][end]
    proportion_IvotedForThisParty = proportionmap(model.properties[:voterBallotTracker][agentid])[myLast_PartyVote]
    proportion_peers_like_me = get_proportion_peers_voteLikeMe(agentid,model)
    # FIXME: check if this proportion is correct

    if HaveIVotedAgainstMyParty(agentid,model)
        keep_party_id_prob = proportion_IvotedForThisParty *  (1-proportion_peers_like_me)
    else
        keep_party_id_prob = proportion_IvotedForThisParty * proportion_peers_like_me
    end

    if (in(agentid,model.properties[:parties_ids])) || (model[agentid].amIaCandidate)
        keep_party_id_prob = 1.0
    end

    return(keep_party_id_prob)
end

function update_partyid!(agentid,model,
                         keep_party_id_prob = get_keep_party_id_prob(agentid,model))
    myLast_PartyVote = model.properties[:voterBallotTracker][agentid][end]

    if rand() > keep_party_id_prob
        model[agentid].myPartyId = myLast_PartyVote
        model.properties[:party_switches][end]+=1
    end
end


function update_streakCounter!(m)
    if m.properties[:incumbent_party] != m.properties[:incumbent_streak_counter].old_incumbentholder
        push!(m.properties[:incumbent_streak_counter].has_switchedlist, true)
        m.properties[:incumbent_streak_counter].current_streak = 1
    else
        m.properties[:incumbent_streak_counter].current_streak +=1
        push!(m.properties[:incumbent_streak_counter].has_switchedlist, false)
    end

    if m.properties[:incumbent_streak_counter].current_streak > m.properties[:incumbent_streak_counter].longest_streak[:streak_value]
        m.properties[:incumbent_streak_counter].longest_streak[:incumbent_pos] = m[m.properties[:incumbent_party]].pos
        m.properties[:incumbent_streak_counter].longest_streak[:streak_value] = m.properties[:incumbent_streak_counter].current_streak
    end
end

# this only makes sense AFTER updating the agent partyid
# but BEFORE updating the global voters_partyids
function add_partyswitch_tocounter!(i,m)
    if m[i].myPartyId != m.properties[:voters_partyids][i]
        m.properties[:party_switches][end]+=1
        end
end

function get_new_supporters(model, old_supporters)
    current_supporters = get_parties_supporters(model)
    new_supporters = kvdictmap((k,v)-> setdiff(v,old_supporters[k]), current_supporters)
    return(new_supporters)
end

function get_mean_among_supporters(supporters, model)
        [mean(model[i].pos[issue]
              for i in supporters)
         for issue in 1:model.properties[:nissues]]
end



function get_new_parties_poss(model, new_supporters, old_supporters)
    ω = model.properties[:ω]

    mean_previous_supporters = dictmap(v->get_mean_among_supporters(v,model),
                                       old_supporters )

    mean_new_supporters = dictmap(v->get_mean_among_supporters(v,model),
                                  new_supporters)

    kvdictmap((k,v)-> (ω .* v .+ ((1-ω) .* mean_new_supporters[k])) |> Tuple, mean_previous_supporters )

end

function set_new_parties_poss!(model,newpposs)
    for (k,v) in newpposs
        model[k].pos = v
        model.properties[:parties_candidateid_ppos_δ][k][:partyposition] = v
    end
end



function set_agent_new_κ!(agentid,model)
    myLast_PartyVote = model.properties[:voterBallotTracker][agentid][end]
    proportion_IvotedForThisParty = proportionmap(model.properties[:voterBallotTracker][agentid])[myLast_PartyVote]
    baseline_κ = model.properties[:κ]
    model[agentid].κ = proportion_IvotedForThisParty * baseline_κ
end

function set_agents_new_κ!(model, kappa_switch= :off)
    for i in abm.allids(model)
        if kappa_switch == :off
            continue
        else
            set_agent_new_κ!(i, model)
        end
    end
end




function model_step!(model)
    candidates_iteration_setup!(model)

    old_supporters = get_parties_supporters(model) |> copy
    set_agents_new_κ!(model, model.properties[:kappa_switch])
    new_winner_party = model[getmostvoted(model, :iteration)].myPartyId

    model.properties[:incumbent_streak_counter].old_incumbentholder = model.properties[:incumbent_party]
    model.properties[:incumbent_party] = new_winner_party
    update_streakCounter!(model)
    for i in abm.allids(model)
        party_i_votedfor = get_whichCandidatePartyAgentVotesfor(i,model)[2]
        push!(model.properties[:voterBallotTracker][i], party_i_votedfor )
    end
    model.properties[:withinpartyshares] = get_withinpartyshares(model)

    push!(model.properties[:party_switches], 0)
    push!(model.properties[:keep_probs], [])
    #=In this loop agents deal with their new choice
    #of candidate by updating their partyid =#
    for i in abm.allids(model)
        keep_prob = get_keep_party_id_prob(i,model)
        push!(model.properties[:keep_probs][end], keep_prob)
        update_partyid!(i,model,keep_prob)
        model.properties[:voters_partyids][i] = model[i].myPartyId

    end

    newposs = get_new_parties_poss(model,get_parties_supporters(model), old_supporters )
    set_new_parties_poss!(model,newposs)
end


# ** Data Collection
#=

FIXME: Test what happens with the  proportion of voters who voted against PartyId candidate

Maybe also some measures of the distribution? Who knows....

• for how long candidate remains winning.
• The proportion of voters voted for someone who was not their party
candidate.
• What is the distance between the candidate the agent voted for and
the candidate of their party.
=#


#adata = [(a->(HaveIVotedAgainstMyParty(a,m)), +)]

function get_distance_MyCandidatevsPartyCandidate(agentid, model)


    closest_to_me = get_closest_candidate(agentid,model)[1]
    mypartycandidate = model.properties[:parties_candidateid_ppos_δ][model[agentid].myPartyId][:partycandidate]
    two_candidates_distance = dist.euclidean(model[closest_to_me].pos,
                                             model[mypartycandidate].pos)
    return(two_candidates_distance)

end

function get_distance_IvsPartyCandidate(agentid,model)
    mypartycandidate = model.properties[:parties_candidateid_ppos_δ][model[agentid].myPartyId][:partycandidate]
    return(dist.euclidean(model[agentid].pos,
                          model[mypartycandidate].pos))
end


get_distance_IvsPartyCandidate(agent::Voter,model) = get_distance_IvsPartyCandidate(agent.id,model)
get_distance_MyCandidatevsPartyCandidate(agent::Voter, model) = get_distance_MyCandidatevsPartyCandidate(agent.id, model)


function get_representativeness(m::abm.ABM)
    -sum([get_distance_IvsPartyCandidate(i,m) for i in abm.allids(m)])/m.properties[:nagents]
end

function get_representativeness(v,m)
    -sum(v)/m.properties[:nagents]
end


function get_partyshare(m)
    proportionmap([m.properties[:voterBallotTracker][agentid][end]
                   for agentid in abm.allids(m)])
end

function get_ENP(m)
    shares = get_partyshare(m)
    return(1/sum([i^2 for i in collect(values(shares))]))
end


function normalized_ENP(m)
    get_ENP(m)/m.properties[:ncandidates]
end


# x->x[x.properties[:incumbent_party]].myPartyId, # get incumbent_party

function get_incumbent_eccentricity(m)
      dist.euclidean(m[m.properties[:incumbent_party]].pos,
                   m.properties[:median_pos])
end

function get_mean_contestant_eccentricity(m)
    contestants = filter(pid-> pid != m.properties[:incumbent_party],
                         m.properties[:parties_ids])
    (contestants .|>
      (contestant -> dist.euclidean(m[contestant].pos,
                                    m.properties[:median_pos])) |>
                                        mean)
end

# function static_preplot!(ax,m)

#     xs,ys = begin
#         poss = [x[:partyposition]
#                 for x in values(m.properties[:parties_candidateid_ppos_δ])]
#         xs = [x[1] for x in poss]
#         ys = [x[2] for x in poss]
#         xs,ys
#     end

#      # be sure that the teacher will be above students
#     obj= GLMakie.scatter!(xs,ys,
#                           marker = :diamond,
#                           markersize = 25,
#                           color = :green)
#     GLMakie.hidedecorations!(ax)
#     GLMakie.translate!(obj, 0, 0, 5)

# end

function saltellidict(varnames::Vector{String}, bounds)
        Dict("names" => varnames,
             "bounds" => PythonCall.PyList(bounds),
             "num_vars" => length(varnames))
    end


function boundsdict_toparamsdf(varnames, bounds;samplesize= 2^8)
    saltelli = pyimport("SALib.sample.saltelli")

    boundsdict = saltellidict(varnames, bounds)

    problem_array = pyconvert(Array,
                              saltelli.sample(boundsdict, samplesize,
                                              calc_second_order = true ))
    foodf = DF.DataFrame()

    for (index,value) in enumerate(boundsdict["names"])
        setproperty!(foodf,Symbol(value), problem_array[:, index])
    end
    return(foodf)
end


datapath = "../../../data"


function run_analysis_onRow(sim_cons,
                            row,
                            rowIdx,
                            whichIter,
                            threadId)
    m = initialize_model(sim_cons[:nagents], sim_cons[:nissues], row.ncandidates,
                             κ=  row.κ, δ =  row.δ)

    adata = [(a->(HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
             (a->(get_distance_IvsPartyCandidate(a,m)), d -> get_representativeness(d,m))]

    mdata = [normalized_ENP,
           x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
             x-> x.properties[:party_switches][end]/x.properties[:nagents],
             get_incumbent_eccentricity]

    ad,md = abm.run!(m, abm.dummystep, model_step!, sim_cons[:niterations] ;
                       adata= adata,
                       mdata=mdata)#, when = collect_steps, when_model= collect_steps)

    alabels = ["step","Prop¬PartyId", "Representativeness"]
    mlabels = ["step", "NENP", "LongestIStreak", "PartySwitches", "Eccentricity"]
    DF.rename!(ad, alabels)
    DF.rename!(md,mlabels)
    df = hcat(ad, md, makeunique=true)
    DF.select!(df, DF.Not(:step_1))

    CSV.write(joinpath(datapath,
                       "row$(rowIdx)_iter$(whichIter)_thread$threadId.csv"),
              df)
    m = nothing
    df = nothing
end


function run_analysis(sim_cons, dm, threadId)

    ProgressMeter.@showprogress 5 "Simulating..."  for iter in 1:10
         for (rowidx,
                                             rowval) in enumerate(eachrow(dm))
            run_analysis_onRow(sim_cons,
                               rowval,
                               rowidx,
                               iter,
                               threadId)
        end
    end
end


get_mean_col_val(df,col,niterations)= StatsBase.mean(df[end-(niterations-1):end,col])

function get_system_measures(df)
    vars =  names(DF.select(df, DF.Not([:step])))

    Dict(Pair(var, get_mean_col_val(df, var, 20)) for var in vars)

end

readdf(dfname) = CSV.read(joinpath(datapath,dfname),
                          DF.DataFrame)

function system_measure_AtRepetions(whichParametization)
    data = readdir(datapath)

    map(get_system_measures∘readdf,
        filter(x-> ("row$(whichParametization)" == split(x, "_")[1]),
               data))
end


function DictionariesToDataFrame(dictlist)
  ret = Dict()                 #Holds dataframe's columns while we build it
  #Get all unique keys from dictlist and make them entries in ret
  for x in unique([y for x in [collect(keys(x)) for x in dictlist] for y in x])
    ret[x] = []
  end
  for row in dictlist          #Loop through each row
    for (key,value) in ret     #Use ret to check all possible keys in row
      if haskey(row,key)       #Is key present in row?
        push!(value, row[key]) #Yes
      else                     #Nope
        push!(value, nothing)  #So add nothing. Keeps columns same length.
      end
    end
  end
  #Fix the data types of the columns
  for (k,v) in ret                             #Consider each column
    row_type = unique([typeof(x) for x in v])  #Get datatypes of each row
    if length(row_type)==1                     #All rows had same datatype
      row_type = row_type[1]                   #Fetch datatype
      ret[k]   = convert(Array{row_type,1}, v) #Convert column to that type
    end
  end
  #DataFrame is ready to go!
  return DF.DataFrame(ret)
end

function get_ParametizationMeasuresMeans(whichparametrization)
    repetitionsvalues= DictionariesToDataFrame(system_measure_AtRepetions(whichparametrization))
    return(DF.mapcols(StatsBase.mean, repetitionsvalues))
end


# TODO: Check performance
# TODO: Turn all dict calls into simple mutating static arrays



end  # this is where the module ends!!!
