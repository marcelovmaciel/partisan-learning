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
import GLMakie
mutable struct Voter{n} <: abm.AbstractAgent
    id::Int
    pos::NTuple{n,Float64}
    amIaCandidate::Bool
    myPartyId::Int
    #= Note that if a voter is a candidate then its ~myPartyId~ should be the
agent's id. Maybe I'll create an dictionary Int=> Symbol to identify the parties
throughout simulation inspection =#
end

const bounds = (0,10)

@kwdef struct ModelParams
    nagents = 3000
    ncandidates = 10
    nissues = 1
    bounds = bounds
    κ = 1. # this influences whether I'll vote against my party or not
    δ = 3. # this influences the radius of candidates a party samples from
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
    actualpartiesposs = Dict(
        Pair(k,
             Dict(:partyposition => v,
                  :partycandidate => -1 ))
        for (k,v) in partiesposs)
    return(actualpartiesposs)
end

function set_candidates!(partiesposs,
                         model::abm.ABM)
    δ = model.properties[:δ]
    getcandidateid(party)= sample(collect(
        abm.nearby_ids(model[party],
                       model,
                       δ,
                       exact = true)))

    candidatepartypairs = []

    for (pid,pvalue) in partiesposs
        # this allows me use it in the initial condition without any problem
        if model.properties[:incumbent] == 0
            #println(getcandidateid(pvalue[:partyposition]),  " ", pid)
            push!(candidatepartypairs,
                  Pair(getcandidateid(pid),pid))
        # this one helps me to jump the incumbent after the initial condition
        elseif pid == model[model.properties[:incumbent]].myPartyId
                continue
        else
            #println(getcandidateid(pvalue[:partyposition]), " ", pid)
            push!(candidatepartypairs,
                  Pair(getcandidateid(pid),
                       pid))
        end
    end

    candidateids =  Dict(candidatepartypairs)
    # println(candidatepartypairs) # BUG: THIS IS WRONG. The candidate id CANNOT be the same as the party for christ sake

    for (candidateid, pid) in candidateids

        model[candidateid].amIaCandidate = true
        model[candidateid].myPartyId = pid
        partiesposs[pid][:partycandidate] = candidateid
    end
end



function set_candidates!(model::abm.ABM)
    set_candidates!(model.properties[:partiesposs],
                    model::abm.ABM)
end



"get_closest_candidate(agentid::Int, model::abm.ABM)"
function get_closest_candidate(agentid::Int,model)
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
        for k in keys(model.properties[:partiesposs]))
end

function get_withinpartyshares(model)
    #= This function is very important and calculates the proportion
    of vote shares WITHIN the parties. That is, to what party
    the supporters of each party have voted for.
    Naturally, it **ONLY MAKES SENSE** if I update
    model.properties[:voterBallotTracker] before using it!!!!=#
    parties_supporters = get_parties_supporters(model)
    withinpartyvotes = Dict(
        Pair(k,
             map(x-> model.properties[:voterBallotTracker][x][end],v))
        for (k,v) in parties_supporters)
   Dict(Pair(k,proportionmap(v)) for (k,v) in withinpartyvotes)
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
    m.properties[:incumbent_streak_counter].old_incumbentholder = m.properties[:incumbent]
    m.properties[:incumbent_streak_counter].current_streak = 1
    m.properties[:incumbent_streak_counter].longest_streak[:streak_value] = 1
    m.properties[:incumbent_streak_counter].longest_streak[:incumbent_pos] = m[m.properties[:incumbent]].pos
end




# TODO: Check if this is indeed correct
"initialize_model(nagents::Int, nissues::Int, ncandidates::Int)"
function initialize_model(nagents::Int, nissues::Int, nparties;
                          κ = 2.,  δ=1.)
    space = abm.ContinuousSpace(ntuple(x -> float(last(bounds)),nissues))

    # postype = typeof(ntuple(x -> 1.,nissues))

    #=
    There are three main auxiliary model collections:
    - What are the parties ids, and their positions and candidates (partiesposs)
    - What are the voters parties ids NOW (voters_partyids)
    - What are the the voters voting track (votersBallotTracker)
    =#
    voters_partyids = Dict{Int64, Int64}()
    voterBallotTracker = Dict{Int64, Vector{Int}}()
    partiesposs = Dict{}()
    withinpartyshares = Dict{}()
    #=
    I am adding that as a model property to later:
    1) use it for synchronous partyid update
    2) plotting (closer to the sugarscape example)
    =#
    properties = Dict(:incumbent => 0,
                      :nagents => nagents,
                      :nissues => nissues,
                      :ncandidates => nparties,
                      :partiesposs => partiesposs,
                      :δ => δ,
                      :κ => κ,
                      :voters_partyids => voters_partyids,
                      :voterBallotTracker=>voterBallotTracker,
                      :withinpartyshares => withinpartyshares,
                      :incumbent_streak_counter => DummyStreakCounter(),
                      :party_switches => [0])
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
        mypartyid = get_closest_candidate(i,model)[2]
        model[i].myPartyId = mypartyid
    end

    model.properties[:incumbent] = new_incumbent
    initialize_incumbent_streak_counter!(model)
    model.properties[:voters_partyids] = Dict((model[x].id => model[x].myPartyId)
                                       for x in abm.allids(model))
    model.properties[:voterBallotTracker] = Dict((k,[v]) for (k,v) in model.properties[:voters_partyids])
    model.properties[:withinpartyshares] = get_withinpartyshares(model)
    return(model)
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
    model[model.properties[:incumbent]].amIaCandidate = true
end


"candidates_iteration_setup!(model::abm.ABM)"
function candidates_iteration_setup!(m::abm.ABM)
    reset_candidates!(m)
    set_candidates!( m)
end


function get_whichCandidatePartyAgentVotesfor(agentid, model)
    κ = model.properties[:κ]

    closest_to_me_id_pid = get_closest_candidate(agentid,model)

    mypartycandidate = model.properties[:partiesposs][model[agentid].myPartyId][:partycandidate]

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



#=
1. Pick the proportion of iterations that I voted for the party I’m voting
for in this iteration.
2. Pick the proportion of people from my party that voted
different from me;
- I change my party id to this other party with a probability equal to tanh(proportion 1 + proportion 2). =#


function update_partyid!(agentid,model)
    myLast_PartyVote = model.properties[:voterBallotTracker][agentid][end]
    proportion_IvotedForThisParty = proportionmap(model.properties[:voterBallotTracker][agentid])[myLast_PartyVote]
    proportion_peersUnlikeMe = (1-get_proportion_peers_voteLikeMe(agentid,model))
    changechange = tanh(proportion_IvotedForThisParty + proportion_peersUnlikeMe)
    if rand() < changechange
        model[agentid].myPartyId = myLast_PartyVote
    end
end


function update_streakCounter!(m)
    if m.properties[:incumbent] != m.properties[:incumbent_streak_counter].old_incumbentholder
        push!(m.properties[:incumbent_streak_counter].has_switchedlist, true)
        m.properties[:incumbent_streak_counter].current_streak = 1
    else
        m.properties[:incumbent_streak_counter].current_streak +=1
        push!(m.properties[:incumbent_streak_counter].has_switchedlist, false)
    end

    if m.properties[:incumbent_streak_counter].current_streak > m.properties[:incumbent_streak_counter].longest_streak[:streak_value]
        m.properties[:incumbent_streak_counter].longest_streak[:incumbent_pos] = m[m.properties[:incumbent]].pos
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




#FIXME: Double-check if I update the model.properties[:partiesposs][:partycandidate]!!!
# I believe it does update with set_candidates! though. Nevertheless, check
function model_step!(model)
    candidates_iteration_setup!(model)

    new_winner = getmostvoted(model, :iteration)

    model.properties[:incumbent_streak_counter].old_incumbentholder = model.properties[:incumbent]
    model.properties[:incumbent] = new_winner
    update_streakCounter!(model)
    for i in abm.allids(model)
        party_i_votedfor = get_whichCandidatePartyAgentVotesfor(i,model)[2]
        push!(model.properties[:voterBallotTracker][i], party_i_votedfor )
    end
    model.properties[:withinpartyshares] = get_withinpartyshares(model)

    push!(model.properties[:party_switches], 0)
    #=In this loop agents deal with their new choice
    #of candidate by updating their partyid =#
    for i in abm.allids(model)
        update_partyid!(i,model)
        add_partyswitch_tocounter!(i,model)
        model.properties[:voters_partyids][i] = model[i].myPartyId
    end
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


function HaveIVotedAgainstMyParty(agentid::Int, model)
    mypartycandidate = model.properties[:partiesposs][model[agentid].myPartyId][:partycandidate]
    mycandidate = get_whichCandidatePartyAgentVotesfor(agentid, model)[1]
    return(mycandidate != mypartycandidate)
end


HaveIVotedAgainstMyParty(agent::Voter, model) = HaveIVotedAgainstMyParty(agent.id,model)

#adata = [(a->(HaveIVotedAgainstMyParty(a,m)), +)]


function get_distance_MyCandidatevsPartyCandidate(agentid, model)


    closest_to_me = get_closest_candidate(agentid,model)[1]
    mypartycandidate = model.properties[:partiesposs][model[agentid].myPartyId][:partycandidate]
    two_candidates_distance = dist.euclidean(model[closest_to_me].pos,
                                             model[mypartycandidate].pos)
    return(two_candidates_distance)

end

function get_distance_IvsPartyCandidate(agentid,model)
    mypartycandidate = model.properties[:partiesposs][model[agentid].myPartyId][:partycandidate]
    dist.euclidean(model[agentid].pos,
                   model[mypartycandidate].pos)
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
# x->x[x.properties[:incumbent]].myPartyId, # get incumbent
mdata = [get_ENP,
         x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
         x-> x.properties[:party_switches][end]]

function static_preplot!(ax,m)

    xs,ys = begin
        poss = [x[:partyposition]
                for x in values(m.properties[:partiesposs])]
        xs = [x[1] for x in poss]
        ys = [x[2] for x in poss]
        xs,ys
    end

     # be sure that the teacher will be above students
    obj= GLMakie.scatter!(xs,ys,
                          marker = :diamond,
                          markersize = 25,
                          color = :green)
    GLMakie.hidedecorations!(ax)
    GLMakie.translate!(obj, 0, 0, 5)

end




# TODO add data collection: incumbent eccentricity! maybe also a mean non-incumbent eccentricity
# TODO: try to visualize how myPartyId evolves!
# TODO: withipartyshares could be an evolving barplot
# TODO: write my own fucking visualization loop

end  # this is where the module ends!!!
