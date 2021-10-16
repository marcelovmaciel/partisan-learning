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

mutable struct Voter{n} <: abm.AbstractAgent
    id::Int
    pos::NTuple{n,Float64}
    amIaCandidate::Bool
    myPartyId::Int
    #= Note that if a voter is a candidate then its ~myPartyId~ should be the
agent's id. Maybe I'll create an dictionary Int=> Symbol to identify the parties
throughout simulation inspection =#
end

const bounds = (0,1)

@kwdef struct ModelParams
    nagents = 3000
    ncandidates = 10
    nissues = 1
    bounds = bounds
    κ = 0.1 # this influences whether I'll vote against my party or not
    ρ = 0.2 # this is the radius of neighbors I take a look. I believe it can be dropped
    ϕ = 0.2 # this is related to neighbors influence (I Think it can be dropped)
    δ = 0.1 # this influences the radius of candidates a party samples from
end


function Voter(id::Int,nissues =  1,
               pos = Tuple(rand(distri.Uniform(bounds...),nissues)))
    amIaCandidate = false
    myPartyId = -3
    return(Voter{nissues}(id,pos,amIaCandidate, myPartyId))
end

function sample_parties_pos(nparties, model)
    ids = collect(abm.allids(model))
    partiesposs = Dict(map(x-> (x,model[x].pos),
                      sample(ids,nparties, replace = false)))
    return(partiesposs)
end


function set_candidates!(partiesposs,model::abm.ABM, δ = 0.1)

    candidateids = Dict(Pair(pid, sample(collect(abm.nearby_ids(pvalue,model,δ, exact = true)))) for (pid,pvalue) in partiesposs)

    for (pid,candidateid) in candidateids
        model[candidateid].amIaCandidate = true
        model[candidateid].myPartyId = pid
        end
end


function set_candidates!(model::abm.ABM)
    set_candidates!(model.properties[:partiesposs],
                    model::abm.ABM,
                    model.properties[:δ])
end



"get_closest_candidate(agentid::Int, model::abm.ABM, mypos_or_partypos::Symbol)"
function get_closest_candidate(agentid::Int,model, mypos_or_partypos = :pos)
    #=
    In this function the optional argument must the field of position that the function must calculate! Either the agents' position or its partyid position!
    =#
    #=
    Those are dummy variables. In the first loop iteration
    they are replaced. Since the bounds: (0,1) there is no
    way distance > dummydistance! the id is negative because
    there can be no such thing. Thus, it will warn me downstream
    if any mistake has happened here.
    =#
    dummydistance = 100.
    dummyid = -2
    for i in abm.allids(model)     #pid.abm.nearby_ids(model[testid].pos, model)
        if model[i].amIaCandidate
            candidatepos = getfield(model[agentid],
                                     mypos_or_partypos)
            distance = dist.euclidean(candidatepos,
                                      model[i].pos)
            if distance < dummydistance
                dummydistance = distance
                dummyid = i
            end
        end
    end
    return(dummyid)
end

"getmostvoted(model::abm.ABM)"
function getmostvoted(model::abm.ABM, initial_or_iteration = :initial)

    # HACK: this ifelse is for reusing this function in the iteration
    # Maybe later write a tagging type for that
    if initial_or_iteration == :initial
        votingfn = get_closest_candidate
    else
        votingfn = get_whoAgentVotesfor
    end

    #=The candidate who is closest to most wins.
    Standard downsian assumption.=#
    closest_candidates = Array{Int}(undef, abm.nagents(model))
    for id in abm.allids(model)
        closest_candidates[id] = votingfn(id,model)
    end
   return(argmax(proportionmap(closest_candidates))) # argmax will return only one maximal, beware of that!
end


# TODO: Check if this is indeed correct
"initialize_model(nagents::Int, nissues::Int, ncandidates::Int)"
function initialize_model(nagents::Int, nissues::Int, nparties,
                          κ = 0.1, ρ = 0.2, ϕ = 0.2, δ=0.1)
    space = abm.ContinuousSpace(ntuple(x -> float(last(bounds)),nissues))

    postype = typeof(ntuple(x -> 1.,nissues))

    #=
    There are three main auxiliary model collections:
    - What are the parties ids and their positions!
    - What are the voters parties ids NOW
    - What are the the voters voting track
    =#


    voters_partyids = Dict{Int64, Int64}()
    voterBallotTracker = Dict{Int64, Vector{Int}}()
    partiesposs = Dict{Int64, postype}()
    #=
    I am adding that as a model property to later:
    1) use it for synchronous partyid update
    2) plotting (closer to the sugarscape example)
    =#
    properties = Dict(:incumbent => 0,
                      :nagents => nagents,
                      :nissues => nissues,
                      :partiesposs => partiesposs,
                      :δ => δ,
                      :κ => κ,
                      :ρ => ρ,
                      :ϕ => ϕ,
                      :voters_partyids => voters_partyids,
                      :voterBallotTracker=>voterBallotTracker)

    model = abm.ABM(Voter{nissues}, space, properties = properties)

    for i in 1:nagents
        vi = Voter(i, nissues)
        abm.add_agent_pos!(vi, model)
    end

    model.properties[:partiesposs] = sample_parties_pos(nparties,
                                                        model)
    set_candidates!(model)
    new_incumbent = getmostvoted(model)

    for i in abm.allids(model)
        mypartyid = begin
            closest_candidate_id = get_closest_candidate(i,model)
            model[closest_candidate_id].myPartyId
        end
        model[i].myPartyId = mypartyid
    end

    model.properties[:incumbent] = new_incumbent
    model.properties[:voters_partyids] = Dict((model[x].id => model[x].myPartyId)
                                       for x in abm.allids(model))
    model.properties[:voterBallotTracker] = Dict((k,[v]) for (k,v) in model.properties[:voters_partyids])
    return(model)
end


# ** Stepping
#=
there is some sublety here.
#basicamente na iteracao vai ter dois tipos de calculo de distancia pra cada
agente ele quer saber quem ele ta mais proximo em termos de posicao DELE MAS
tambem quem que ta mais proximo da partyid dele entao ele considera dois
candidatos se o candidato mais proximo dele 'e tipo muito mais proximo que o
candidato do partido, uma contante kappa, ele vota no outro maas, ao faze-lo sua
partyid passa a ser a media da posicao desse novo candidato com a partyid
anterior
=#
"reset_candidates!(model::abm.ABM)"
function reset_candidates!(model::abm.ABM)
    for agent in abm.allids(model)
        model[agent].amIaCandidate = false
    end
    model[model.properties[:incumbent]].amIaCandidate = true
end

"candidates_iteration_setup!(model::abm.ABM)"
function candidates_iteration_setup!(m::abm.ABM)
    reset_candidates!(m)
    set_candidates!(m.properties[:ncandidates]-1, m)
end

"get_whoAgentVotesfor(agentid, model)"
function get_whoAgentVotesfor(agentid, model)
    κ = model.properties[:κ]
    closest_to_me = get_closest_candidate(agentid,model)
    closest_to_myPartyId = get_closest_candidate(agentid,
                                                 model,
                                                 :myPartyId)
    whoillvotefor = closest_to_myPartyId
    two_candidates_distance = dist.euclidean(model[closest_to_me].pos,
                                             model[closest_to_myPartyId].pos)

    if two_candidates_distance > κ
        whoillvotefor = closest_to_me
    end

    return(whoillvotefor)

end


function update_partyid!(agentid,model)
    # let modelsize = 1000
    whoIllVote_now = get_whoAgentVotesfor(agentid, model) #this guy allocates  17 times!
    model[agentid].myPartyId .= broadcast((x,y) -> mean([x,y]), model[agentid].myPartyId,
                                                        model[whoIllVote_now].pos) # this guy allocates 3 times!
    model.properties[:partyids][agentid] = model[agentid].myPartyId
end


"function get_mean_neighborhoodPid(agentid,model, ρ)
This is a helper for updatePid_neighbors_influence!"
function get_mean_neighborhoodPid(agentid,model, ρ)
    neighbors = abm.nearby_ids(model[agentid].pos, model, ρ )
    #mean_neighbor_Pid = MVector{n}(ntuple(x-> 0., n)), do this later
    StatsBase.mean(map(x->model.properties[:partyids][x], neighbors))
    #= Notice that I'm looking at the dictionary!
    With that I can mutate synchronously the voters myPartyIds.
    After I mutate all of them I gotta udpate the model.properties[:partyids]=#
    #= FIXME: This map allocates a lot. Later I'll fix that by rolling my mean loop =#
end



function updatePid_neighbors_influence!(agentid,model, ρ, ϕ)
    if rand() < ϕ
        neighborsmean = get_mean_neighborhoodPid(agentid,model, ρ)
        model[agentid].myPartyId .= broadcast((x,y) -> mean([x,y]),
                                              model[agentid].myPartyId,
                                              neighborsmean)
    end
end


function model_step!(model)
    candidates_iteration_setup!(model)

    new_winner = getmostvoted(model, :iteration)

    model.properties[:incumbent]= new_winner

    #=In this loop agents deal with their new choice
    #of candidate by updating their partyid =#
    for i in abm.allids(model)
        update_partyid!(i,model)
    end

    # In this loop agents deal with their neighbors' udpates
    for i in abm.allids(model)
        updatePid_neighbors_influence!(i,model,
                                       model.properties[:ρ],
                                       model.properties[:ϕ])
    end

    for i in abm.allids(model)
        model.properties[:partyids][i] = model[i].myPartyId
    end

end

# ** Data Collection

#=

FIXME: Test what happens with the  proportion of voters who voted against PartyId candidate
# Maybe also some measures of the distribution? Who knows....

=#

function HaveIVotedAgainstMyParty(agentid::Int, model)
    closest_to_myPartyId = get_closest_candidate(agentid,
                                                 model,
                                                 :myPartyId)
    return(get_whoAgentVotesfor(agentid, model) != closest_to_myPartyId)
end

HaveIVotedAgainstMyParty(agent::Voter, model) = HaveIVotedAgainstMyParty(agent.id,model)

#adata = [(a->(HaveIVotedAgainstMyParty(a,m)), +)]
mdata = [:incumbent]


function get_distance_IvsParty(agentid, model)
    κ = model.properties[:κ]
    closest_to_me = get_closest_candidate(agentid,model)
    closest_to_myPartyId = get_closest_candidate(agentid,
                                                 model,
                                                 :myPartyId)
    whoillvotefor = closest_to_myPartyId
    two_candidates_distance = dist.euclidean(model[closest_to_me].pos,
                                             model[closest_to_myPartyId].pos)
    return(two_candidates_distance)

end

get_distance_IvsParty(agent::Voter, model) = get_distance_IvsParty(agent.id, model)

# TODO add data collection: incumbent eccentricity! maybe also a mean non-incumbent eccentricity


end  # this is where the module ends!!!
