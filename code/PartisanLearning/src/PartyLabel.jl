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
# ** Initial condition

const bounds = (0,100)

mutable struct Voter{n} <: abm.AbstractAgent
    id::Int
    pos::NTuple{n,Float64}
    amIaCandidate::Bool
    κ::Float64
    myPartyId::Int
end

function sample_uniform_pos(nissues = 2)
Tuple(rand(distri.Uniform(bounds...),nissues))
end

OVL(test_cohend) = 2*distri.cdf(distri.Normal(),
                                -abs(test_cohend)/2)
one_modal_dispersed = (distri.Normal(50,25),
                       distri.Normal(50,25))

overlap_50_poss = (distri.Normal(43.25,10), distri.Normal(56.75,10))
overlap_80_poss = (distri.Normal(47.5,10),distri.Normal(52.5,10) )
overlap_20_poss = (distri.Normal(37.25, 10), distri.Normal(62.75,10))


more_dispersed_1d_poss = (distri.Normal(50 - 20.25/2,15),
                          distri.Normal(50 + 20.25/2,15))

function sample_1dnormal(bound, poss = overlap_50_poss)

    #= - I'll simply put a constant in one dimension and work on the other
    lol.ap - I'll use cohen d to define overlapping distributions: - As shown
    here [[https://rpsychologist.com/cohend/][Interpreting Cohen&#x27;s d | R
    Psychologist]] a cohen d of 1.35 gives an overlap of 50%. - Thus if I have a
    distribution of (43,10) I gotta have the other as (56.5,10)
    δ= μ2​−μ1/σ​​ is the formula for a cohend.
    =#

    if rand([false,true])
        pos = (rand(poss[1]), bound[2]/2 )
        if pos[1] < 0.0
            pos = (0.01, bound[2]/2 )
        elseif pos[1] > 100
            pos = (99.999, bound[2]/2 )
        end

    else
        pos = (rand(poss[2]),bound[2]/2 )
        if pos[1] > 100
            pos = (99.999, bound[2]/2 )
        elseif pos[1] < 0.0
            pos = (0.01, bound[2]/2 )
            end
    end
    return(pos)
end


overlap_initializor(whichdist) = () -> sample_1dnormal((100., 5.), whichdist)


function Voter(id::Int;nissues =  1,
               pos = () -> sample_uniform_pos(nssissues), κ = 10.)
    amIaCandidate = false
    myPartyId = -3
    return(Voter{nissues}(id,pos(),amIaCandidate, κ, myPartyId))
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

function simple_sample_candidates(party, m, n=1)
    sample(
        collect(abm.nearby_ids(m[party], m, m.properties[:δ],
                                  exact = true)),n)
end

function simple_select_primariesCandidates(model::abm.ABM, n =1)

    party_candidate_pairs = (
        model.properties[:parties_ids] .|>
            (pid -> Pair(pid,simple_sample_candidates(pid,model, n))) |>
            Dict)

      return(party_candidate_pairs)
end


function sample_candidates2(party,m,n=4)
    δ = m.properties[:δ]

    nearby_agents = collect(abm.nearby_ids(m[party], m, δ, exact = true))

    if m.properties[:is_at_step] >=4


    nearby_agents_from_the_party = filter(agentid -> m[agentid].myPartyId == party,
                                          nearby_agents)

    turnout_supporters = get_supporters_who_turnout(nearby_agents_from_the_party, m)

        (isempty(turnout_supporters) ?
        fill(party, (n,1)) :
        sample(turnout_supporters,n))

    else
    nearby_agents_from_the_party = filter(agentid -> m[agentid].myPartyId == party,
                                          nearby_agents)
    (isempty(nearby_agents_from_the_party) ?
        fill(party, (n,1)) :
        sample(nearby_agents_from_the_party,n))
    end

end

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

function sample_candidates(party, m)
    if m.properties[:is_at_step] >=4

        party_supporters = get_parties_supporters(m)[party]

        turnout_supporters = get_supporters_who_turnout(party_supporters,
                                                    m)
        (isempty(turnout_supporters) ?
            fill(party, (4,1)) :
            sample(turnout_supporters,4))

    else

        nearby_agents = collect(abm.nearby_ids(m[party],
                                               m, m.properties[:δ],
                                               exact = true))
        (isempty(nearby_agents) ?
            fill(party, (4,1)) :
            sample(nearby_agents,4))
    end

end


function select_primariesCandidates(model::abm.ABM)
    party_candidate_pairs = (
        model.properties[:parties_ids] .|>
            (pid -> Pair(pid,sample_candidates(pid,model))) |>
            Dict)

      return(party_candidate_pairs)
end

function get_plurality_result(primariesresult::Dict)

    pm =  dictmap(proportionmap, primariesresult)
    if any(x-> length(x) == 0, values(pm))
        #println(primariesresult, "defuck is here")
        newresult = Dict((k, length(v) == 0 ? [k] : v ) for (k,v) in primariesresult)

           #println(newresult)
           pm = dictmap(proportionmap,newresult)
    end
    d = dictmap(argmax, pm)
    return(d)
end

function get_plurality_result(m::abm.ABM, seconditer_switch=false)
    if seconditer_switch
        primariesresult = secondIt_get_primaries_votes(m)
        result = get_plurality_result(primariesresult)
    else

        candidates = select_primariesCandidates(m)
        primariesresult = get_primaries_votes(m,candidates)


        pm = dictmap(proportionmap, primariesresult)

       if any(x-> length(x) == 0, values(pm))
           #println(primariesresult, "defuck is here2", pm)

           newresult = Dict((k, length(v) == 0 ? [k] : v ) for (k,v) in primariesresult)

           #println(newresult)
           pm = dictmap(proportionmap,newresult)
           println(map(x-> length(x), values(pm)))

       end
        result = dictmap(argmax, pm)
        #result = get_plurality_result(primariesresult)
    end
    return(result)
end

# FIXME: this is also wrong. I'm so tired of this code...
function get_runoff_result(;primariesresult=primariesresult,m=m, seconditer_switch = false)

    #println(m)
    if seconditer_switch
        #println(primariesresult)
        primariesproportion = dictmap(proportionmap
                                      ,primariesresult)
        function runoff(k,v)
            if any(x->x>0.5,values(v))
            #println(v)
            argmax(v)
            else
             #    println(secondIt_get_parties_supporters(m)[k])
                 toptwo = sort(collect(v),
                               by=x->x[2],
                               rev = true)[1:2] .|> x->x[1]
                 (map(x->get_closest_fromList(x,toptwo,m),
                  secondIt_get_parties_supporters(m)[k]) |>
                      proportionmap |>
                      argmax)
            end
        end

        return(kvdictmap(runoff, primariesproportion))

    else

        primariesproportion = dictmap(proportionmap
                                      ,primariesresult)

    function runoff2(k,v)
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
        return(kvdictmap(runoff2, primariesproportion))
    end
end

function get_runoff_result(m, seconditer_switch = false)
    if seconditer_switch

    primariesresult = secondIt_get_primaries_votes(m)


    runoffresult = get_runoff_result(primariesresult=primariesresult,m=m, seconditer_switch=seconditer_switch)

    else
        candidates = select_primariesCandidates(m)

        primariesresult = get_primaries_votes(m,candidates)
        #println(primariesresult)
        runoffresult =  get_runoff_result(primariesresult=primariesresult,m=m, seconditer_switch=seconditer_switch)

    end

    return(runoffresult)

end


function set_candidates!(model, switch)

    if switch == :initial

        candidateids = simple_select_primariesCandidates(model)

    elseif switch == :second
        if model.properties[:switch] == :runoff
            candidateids = get_runoff_result(model,
                                             true)
        elseif model.properties[:switch] == :plurality
            candidateids = get_plurality_result(model, true)
        end
    elseif switch == :runoff
        candidateids = get_runoff_result(model,
                                        false)
    elseif switch == :plurality
        candidateids = get_plurality_result(model)
    else
        candidateids = get_random_candidates(model)
    end


    for (pid, candidateid) in candidateids
        if typeof(candidateid) <: Array
            candidateid = first(candidateid)
        end
        # println(candidateids,"in set candidates" )
        model[candidateid].amIaCandidate = true

        model[candidateid].myPartyId = pid

        # TODO: i forgot to update that somewhere!
        model.properties[:parties_candidateid_ppos_δ][pid][:partycandidate] = candidateid

    end
#println(candidateids,"in set candidates2" )
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

function secondIt_get_parties_supporters(model)
    overall_placeholder = []
    for k in model.properties[:parties_ids]
        placeholder = []
        for agent in abm.allids(model)
            if model.properties[:voterBallotTracker][agent][end] == k
                push!(placeholder, agent)
            end
            push!(overall_placeholder, Pair(k, placeholder))
        end
    end
        return(Dict(overall_placeholder))
end


function will_I_turnout(i,m)

    if >(m.properties[:is_at_step],1) &&  <(m.properties[:is_at_step], 4)
        lastvote = m.properties[:voterBallotTracker][i][end]
        will_I = <(rand(distri.Uniform(0.75,1)),
                   proportionmap(m.properties[:voterBallotTracker][i])[lastvote])

    else

        voter = m[i]

        voterballot = m.properties[:voterBallotTracker][i]

        proportion_votesi = proportionmap(voterballot)


        if !(voter.myPartyId in
         keys(proportion_votesi))
        proportionvoted_for_party = 0.0
        else
            proportionvoted_for_party = proportion_votesi[voter.myPartyId]
        end

        will_I = rand(distri.Uniform(0.75,1)) < proportionvoted_for_party
    end

    return(will_I)
end


function get_supporters_who_turnout(supporters,m)
    filter(i->will_I_turnout(i,m),supporters)
end


function get_parties_supporters(model)
    Dict(
        Pair(k,
             collect(keys(filter(t->t[2]==k,
                    model.properties[:voters_partyids]))))
        for k in model.properties[:parties_ids])
end


function secondIt_get_primaries_votes(m,
                                      primariesCandidatesDict = simple_select_primariesCandidates(m,4))

    all_parties_supporters = secondIt_get_parties_supporters(m)
    parties_supporters = dictmap(supporters -> get_supporters_who_turnout(supporters, m),
                                 all_parties_supporters)
    get_closest_toI(i) = get_closest_fromList(i,
                         primariesCandidatesDict[m.properties[:voterBallotTracker][i][end]],
                                              m)

    return(dictmap(supporters->map(get_closest_toI,supporters),
                parties_supporters))
end

function get_primaries_votes(m, primariesCandidatesDict)

    all_parties_supporters = get_parties_supporters(m)

    parties_supporters = dictmap(supporters -> get_supporters_who_turnout(supporters, m),
                                 all_parties_supporters)

    # for (k,v) in parties_supporters
    #     if length(v) == 0
    #         parties_supporters[k] = [k]
    #     end
    # end

    # FIXME: im getting an error here, ghe length is all agents !!!
    # it shouldn't be
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

@kwdef struct ModelParams
    nagents = 3000
    nparties = 2
    nissues = 2
    bounds = bounds
    κ = 1. # this influences whether I'll vote against my party or not
    δ = 3. # this influences the radius of candidates a party samples from
    ω = 0.8
    switch = :plurality
    kappa_switch = :off
    special_bounds = (true, (100., 5.))
    party_pos_hardwired = false
    voter_pos_initializor = overlap_initializor(overlap_50_poss)
end


function initialize_model(;nagents = 100,
                          nissues = 2,
                          nparties = 2,
                          κ = 2.,
                          δ=1.,
                          switch= :random,
                          bounds = bounds,
                          ω = 0.8,
                          kappa_switch = :off,
                          special_bounds = (false, bounds),
                          voter_pos_initializor = () -> sample_1dnormal(special_bounds[2]),
                          party_pos_hardwired = false)

    if special_bounds[1] == true
        space = abm.ContinuousSpace(special_bounds[2], periodic = false)
    else
        space = abm.ContinuousSpace(ntuple(x -> float(last(bounds)),nissues), periodic = false)
    end


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
                      :cross_voting => [0],
                      :median_pos => [],
                      :kappa_switch => kappa_switch,
                      :is_at_step => 0)

    model = abm.ABM(Voter{nissues}, space;properties = properties)

    for i in 1:nagents
        vi = Voter(i, nissues=nissues, κ = model.properties[:κ], pos = voter_pos_initializor)
        abm.add_agent_pos!(vi, model)
    end

    if !party_pos_hardwired
    model.properties[:parties_candidateid_ppos_δ] = sample_parties_pos(nparties,
                                                                       model)
    else
        model[1].pos = (15, 5/2)
        model[2].pos = (85, 5/2)

        partiesposs = Dict(map(x-> (x,model[x].pos),[1,2]))

        actualpartiesposs = dictmap(v-> Dict(:partyposition => v,
                                         :partycandidate => -1,
                                         :δ => model.properties[:δ]),
                                    partiesposs)

        model.properties[:parties_candidateid_ppos_δ] = actualpartiesposs

    end

    for (i,v) in enumerate(collect(keys(model.properties[:parties_candidateid_ppos_δ])))
        model.properties[:parties_ids][i]=v
    end

    # const model.properties[:parties_ids] = Tuple(model.properties[:parties_ids])

    for id in model.properties[:parties_ids]
        model[id].myPartyId = id
    end


    # for i in abm.allids(model)
    #     mypartyid = get_closest_fromList(i,model.properties[:parties_ids],model)
    #     model[i].myPartyId = mypartyid
    # end

    model.properties[:voters_partyids] = Dict(
        (model[x].id => -1)
        for x in abm.allids(model))

    initialize_incumbent_streak_counter!(model)
    model.properties[:voterBallotTracker] = Dict((k,Int64[]) for (k,_) in model.properties[:voters_partyids])
    model.properties[:median_pos] = get_median_pos(model)
    return(model)
end


initialize_model(x::ModelParams) = initialize_model(;ntfromstruct(x)...)

function assume_initial_partyid!(i, m)
    if (m.properties[:is_at_step] == 4) && (!in(i,m.properties[:parties_ids]))
    pairs = collect(proportionmap(m.properties[:voterBallotTracker][i]))
    vals = first.(pairs)
    weights = last.(pairs)
    which_party_id_to_assume = sample(vals, Weights(weights))
    m[i].myPartyId = which_party_id_to_assume
    m.properties[:voters_partyids][i] = m[i].myPartyId
    else
        nothing
    end
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
function candidates_iteration_setup!(m::abm.ABM, switch)
    reset_candidates!(m)
    set_candidates!(m, switch)
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


function add_crossvoting_tocounter!(i,m)
    if HaveIVotedAgainstMyParty(i,m)
        m.properties[:cross_voting][end] +=1
    end
end


function get_new_supporters(model, old_supporters)
    current_supporters = get_parties_supporters(model)
    new_supporters = kvdictmap((k,v)-> setdiff(v,old_supporters[k]), current_supporters)
    return(new_supporters)
end

function get_mean_among_supporters(supporters::Vector, model)

        [mean(model[i].pos[issue]
              for i in supporters)
         for issue in 1:model.properties[:nissues]]
end

function get_new_parties_poss(model, new_supporters, old_supporters)
    ω = model.properties[:ω]

    mean_previous_supporters = dictmap(v->get_mean_among_supporters(v,model),
                                       old_supporters )
    # new_supporters = kvdictmap((k,v)-> isempty(v) ? old_supporters[k] : v, new_supporters)
    # if any(values(dictmap(isempty,new_supporters)))
    #               mean_new_supporters = mean_previous_supporters
    # else
    mean_new_supporters = kvdictmap((k,v)-> isempty(v) ? model[k].pos : get_mean_among_supporters(v,model),
                                  new_supporters)
    #end
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

    model.properties[:is_at_step] += 1
    if model.properties[:is_at_step] == 1
        candidates_iteration_setup!(model, :initial)
        for i in abm.allids(model)
            candidates = map(x->x[:partycandidate], values(model.properties[:parties_candidateid_ppos_δ]))

            myvote = get_closest_fromList(i,candidates,model)
            # FIXME: something broken here
            push!(model.properties[:voterBallotTracker][i], model[myvote].myPartyId)

        end

    elseif (model.properties[:is_at_step] == 2) || (model.properties[:is_at_step] == 3)
        candidates_iteration_setup!(model, :second)
        for i in abm.allids(model)
            candidates = map(x->x[:partycandidate], values(model.properties[:parties_candidateid_ppos_δ]))
            myvote = get_closest_fromList(i,candidates,model)
            push!(model.properties[:voterBallotTracker][i], model[myvote].myPartyId)
        end

    else

        for i in abm.allids(model) assume_initial_partyid!(i, model) end
        candidates_iteration_setup!(model, model.properties[:switch])

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
        push!(model.properties[:cross_voting], 0)

        for i in abm.allids(model)
            add_crossvoting_tocounter!(i, model)
        end

        push!(model.properties[:keep_probs], [])
        #=In this loop agents deal with their new choice
        #of candidate by updating their partyid =#
        for i in abm.allids(model)
            keep_prob = get_keep_party_id_prob(i,model)
            push!(model.properties[:keep_probs][end], keep_prob)
            update_partyid!(i,model,keep_prob)
            model.properties[:voters_partyids][i] = model[i].myPartyId

        end

        newposs = get_new_parties_poss(model,
                                       get_new_supporters(model,old_supporters),
                                       old_supporters )

        set_new_parties_poss!(model,newposs)

    end
end

foursteps!(m) =  for _ in 1:4 model_step!(m) end
