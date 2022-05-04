mutable struct Voter{n} <: abm.AbstractAgent
    id::Int
    pos::NTuple{n,Float64}
    amIaCandidate::Bool
    myPartyId::Int
    κ::Float64
end

const bounds = (0, 100)
const oneDBound = (100., 5.)

OVL(test_cohend) = 2*distri.cdf(distri.Normal(),
                                -abs(test_cohend)/2)

const one_modal_dispersed = (distri.Normal(50,25),
                       distri.Normal(50,25))

const overlap_50_poss = (distri.Normal(43.25,10), distri.Normal(56.75,10))
const overlap_80_poss = (distri.Normal(47.5,10),distri.Normal(52.5,10) )
const overlap_20_poss = (distri.Normal(37.25, 10), distri.Normal(62.75,10))

const more_dispersed_1d_poss = (distri.Normal(50 - 20.25/2,15),
                                distri.Normal(50 + 20.25/2,15))
#= 
There are two types of sampling here: 
- Simple global distribution; 
- Overlapping normal dists
=#

function general_sampling(nissues,
                          distribution)
Tuple(rand(distribution,nissues))
end

function sample_uniform_pos(nissues)
    general_sampling(nissues, distri.Uniform(bounds...))
end

function sample_simple_normal(nissues)
    general_sampling(nissues,
    one_modal_dispersed[1])
end

# Note that in this formulation the agent will have npositions all from the
# same dist! 
function sample_overlapping_normal(nissues, poss)
    if rand([false,true])
        d = poss[1] 

    else
        d = poss[2]
    end

    pos = general_sampling(nissues,d)
    
    function correctbounds(p)    
        if p < 0.0
            p = 0.01
        elseif p > 100
            p = 99.999 
        end
        return(p)    
    end    
    pos = correctbounds.(pos)
    return(pos)
end

sample_overlapping_normal(d) = sample_overlapping(1,d)

function sample_overlapping_2d_to_1d_hack(d)
    (first(sample_overlapping_normal(1,d)), oneDBound[2]/2)
end

function Voter(id::Int; nissues =  1,
               pos, κ = 10.)
    amIaCandidate = false
    myPartyId = -3

    v = Voter{nissues}(id, pos, 
       amIaCandidate,    myPartyId,    κ)
    return(v)
end


# ** Sampling parties positions
#=
There are, many ways of sampling parties positions: 
- Simply picking at random; 
- Bias it towards either center or extremes; 
- Hardwire it!  
=#

# this works for any nymber of dimensions
function sample_parties_pos(nparties::Int, model)
    ids = collect(abm.allids(model))

        partiesposs = Dict(map(x-> (x,model[x].pos),
                           sample(ids,nparties, replace = false)))
    actualpartiesposs = dictmap(v-> Dict(:partyposition => v,
                                         :partycandidate => -1,
                                         :δ => model.properties[:δ]),
                                partiesposs)
    return(actualpartiesposs)
end

function set_parties_hardwired_poss!(hardwired_positions,model)
    for (i,v) in enumerate(hardwired_positions)
        model[i] = v
    end
    harwired_ids = collect(1:length(hardwired_positions))
    pposs = (map(x-> (x, model[x].pos), harwired_ids) |>
        Dict |>
        foo ->  dictmap(v-> Dict(:partyposition => v,
                              :partycandidate => -1,
                                 :δ => model.properties[:δ]), foo))
    return(pposs)
end

# ** Getting Parties Supporters
#= Here there are three ways of getting supporters
- The normal way;
- The pre-step way; 
- The filtering for turnout way; 
=#

# this seems to work, but is ugly as fuck
function get_parties_supporters(model)
    Dict(Pair(k,collect(keys(filter(t->t[2]==k,
                    model.properties[:voters_partyids]))))
        for k in model.properties[:parties_ids])
end


# I don't fucking rememember why I wrote this shit 
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


# FIXME: this should be after the model intitialized!
function will_I_turnout(i,m)

    # how partisan must I be to turnout
    turnout_region = distri.Uniform(0.75,1)

    if >(m.properties[:is_at_step],1) &&  <(m.properties[:is_at_step], 4)
        lastvote = m.properties[:voterBallotTracker][i][end]
        will_I = <(rand(turnout_region),
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

        will_I = rand(turnout_region) < proportionvoted_for_party
    end

    return(will_I)
end


function get_supporters_who_turnout(supporters,m)
    filter(i->will_I_turnout(i,m),supporters)
end

# * Sampling candidates for primaries 


function get_party_randomCandidate_dict(m)
    parties_supporters = get_parties_supporters(m)
        if m.properties[:incumbent_party] == 0
            dictmap(rand,parties_supporters)
    else
        delete!(parties_supporters, m[m.properties[:incumbent_party]].myPartyId)
        dictmap(rand,
                parties_supporters)
    end
end

function initial_steps_candidate_sampling(party,m, n=1)
    nearby_agents = collect(abm.nearby_ids(m[party],m,
     m.properties[:δ], exact = true))

    (isempty(nearby_agents) ? 
    fill(party, (n,1)) :
    sample(nearby_agents,n))

end 


function normal_steps_candidate_sampling(party,m)
    party_supporters = get_parties_supporters(m)[party]

    turnout_supporters = get_supporters_who_turnout(party_supporters,
                                                m)
    (isempty(turnout_supporters) ?
        fill(party, (4,1)) :
        sample(turnout_supporters,4))
end

function sample_candidates(party, m)
    if m.properties[:is_at_step] >=4
        normal_steps_candidate_sampling(party,m)
    else
        initial_steps_candidate_sampling(party,m )
    end
end

function party_candidate_Dict_generator(sampler, model)
    (model.properties[:parties_ids] .|>
    (pid -> Pair(pid,sampler(pid,model))) |>
    Dict)
end

function simple_select_primariesCandidates(model, n) 
    party_candidate_Dict_generator((p,m)->  initial_steps_candidate_sampling(p,m,n), 
    model)
end

function select_primariesCandidates(model::abm.ABM)
  party_candidate_Dict_generator(sample_candidates, model)
end


# * Primaries Election 

# ** Primaries Voting Procedures
# FIXME: This is branchy as fuck. FIX THAT
function get_plurality_result(primariesresult::Dict)

    pm =  dictmap(proportionmap, primariesresult)
    if any(x-> length(x) == 0, values(pm))

        newresult = Dict((k, length(v) == 0 ? [k] : v ) for (k,v) in primariesresult)
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
           newresult = Dict((k, length(v) == 0 ? [k] : v ) for (k,v) in primariesresult)
           pm = dictmap(proportionmap,newresult)
       end
        result = dictmap(argmax, pm)
    end
    return(result)
end


function second_round_runoff_calculus(k,v, m,  get_supporters_fn)
    toptwo = sort(collect(v), by=x->x[2], rev = true)[1:2] .|> x->x[1]
    (map(x->get_closest_fromList(x,toptwo,m),
     get_supporters_fn(m)[k]) |>
     proportionmap |>
     argmax)
end


function seconditer_runoff(primariesresult, m)
   primariesproportion = dictmap(proportionmap
                                      ,primariesresult)
        function runoff(k,v)
            if any(x->x>0.5,values(v))
            argmax(v)
            else
                second_round_runoff_calculus(k,v, m,
                 secondIt_get_parties_supporters)
            end
        end
        return(kvdictmap(runoff, primariesproportion))
end

function normal_iteration_runoff(primariesresult, m)
    primariesproportion = dictmap(proportionmap, primariesresult)
    function runoff(k,v)
        if any(x->x>0.5,values(v))
                argmax(v)
            else
                second_round_runoff_calculus(k,v, m,
                get_parties_supporters)
            end
        end
return(kvdictmap(runoff, primariesproportion))
end


function get_runoff_result(;primariesresult=primariesresult,m=m, seconditer_switch = false)
    if seconditer_switch
        seconditer_runoff(primariesresult, m)

    else
        normal_iteration_runoff(primariesresult, m) 
    end
end


function get_runoff_result(m, seconditer_switch = false)
    if seconditer_switch

    primariesresult = secondIt_get_primaries_votes(m)


    runoffresult = get_runoff_result(primariesresult=primariesresult, m=m, 
    seconditer_switch=seconditer_switch)

    else
        candidates = select_primariesCandidates(m)

        primariesresult = get_primaries_votes(m,candidates)
        runoffresult =  get_runoff_result(primariesresult=primariesresult,m=m,
         seconditer_switch=seconditer_switch)

    end

    return(runoffresult)

end




# FIXME: This is impossible to extend and gotta change
# Val types might help me!  
function set_candidates!(model, switch)

    if switch == :initial

        candidateids = simple_select_primariesCandidates(model, 1)

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
        candidateids = get_party_randomCandidate_dict(model)
    end


    for (pid, candidateid) in candidateids
        if typeof(candidateid) <: Array
            candidateid = first(candidateid)
        end
        model[candidateid].amIaCandidate = true

        model[candidateid].myPartyId = pid

        model.properties[:parties_candidateid_ppos_δ][pid][:partycandidate] = candidateid

    end
end

# ** Getting closest candidate for each voter

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

# ** Getting primaries votes

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
    # it shouldn't be
    get_closest_toI(i) = get_closest_fromList(i,
                         primariesCandidatesDict[m[i].myPartyId],
                                              m)
    return(dictmap(supporters->map(get_closest_toI,supporters),
                parties_supporters))
end

# * Getting who would win

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

# FIXME: this is an within model measure
# TODO: think how to generalize this to any model!

DummyStreakCounter() = StreakCounter(1,
                                     Vector{Bool}(),
                                     0,
                                     Dict(:streak_value => 0,
                                          :incumbent_pos => (0,)))

function initialize_incumbent_streak_counter!(m)
    m.properties[:incumbent_streak_counter].old_incumbentholder = 0
    m.properties[:incumbent_streak_counter].current_streak = 0
    m.properties[:incumbent_streak_counter].longest_streak[:streak_value] = 0
    m.properties[:incumbent_streak_counter].longest_streak[:incumbent_pos] = ntuple(x->0.,Val(m.properties[:nissues]))
end

hold(f,ps) = () -> f(ps)

# * Initializing the model

# TODO: every model type should have one such constructor!
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
    voter_pos_initializor = hold(sample_overlapping_2d_to_1d_hack,
     overlap_50_poss)
end


function get_space_given_special_bound_switch(bounds, special_bounds, nissues)
    if special_bounds[1] == true
        space = abm.ContinuousSpace(special_bounds[2], periodic = false)
    else
        space = abm.ContinuousSpace(ntuple(x -> float(last(bounds)),nissues), periodic = false)
    end
    return(space)
end 


function initialize_m_properties(nparties, nagents, nissues, δ, κ, ω, kappa_switch, switch ) 
    voters_partyids = Dict{Int64, Int64}()
    voterBallotTracker = Dict{Int64, Vector{Int}}()
    parties_candidateid_ppos = Dict{}()
    parties_ids = Vector{Int}(undef,nparties)
    withinpartyshares = Dict{}()

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
    return(properties)
end


function add_agents!(model, nagents, nissues, voter_pos_initializor)
    
    for i in 1:nagents
    

        pos = voter_pos_initializor()
        vi = Voter(i, nissues=nissues, pos = pos,  κ = model.properties[:κ])
        
        abm.add_agent_pos!(vi, model)
    end
end

function set_hardwired_pposs(model)
    model[1].pos = (15, 5/2)
    model[2].pos = (85, 5/2)

    partiesposs = Dict(map(x-> (x,model[x].pos),[1,2]))

    actualpartiesposs = dictmap(v-> Dict(:partyposition => v,
                                     :partycandidate => -1,
                                     :δ => model.properties[:δ]),
                                partiesposs)

    model.properties[:parties_candidateid_ppos_δ] = actualpartiesposs
end

function set_sampled_pposs!(model)
    model.properties[:parties_candidateid_ppos_δ] = sample_parties_pos(model.properties[:ncandidates],
    model) # TODO: fix this function call. Only model is needed     
end    

function set_partyposs(model, party_pos_hardwired)
    if !party_pos_hardwired
        set_sampled_pposs!(model)
    else
        set_hardwired_pposs!(model)
        end
end

function update_pposs_tracker!(model)
    for (i,v) in enumerate(collect(keys(model.properties[:parties_candidateid_ppos_δ])))
        model.properties[:parties_ids][i]=v
    end
end

function update_voters_partyids!(model) # FIXME: is that even necessary?
    for id in model.properties[:parties_ids]
        model[id].myPartyId = id
    end
end

function initialize_vpids_tracker(model)
    model.properties[:voters_partyids] = Dict(
        (model[x].id => -1)
        for x in abm.allids(model))
end


function initialize_ballot_trackers(model)
    model.properties[:voterBallotTracker] = Dict((k,Int64[]) for (k,_) in model.properties[:voters_partyids])
end 

function initialize_median_pos_tracker(model)
    model.properties[:median_pos] = get_median_pos(model)
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

    space = get_space_given_special_bound_switch(bounds, special_bounds, nissues)

    properties = initialize_m_properties(nparties, nagents, nissues, δ, κ, ω, kappa_switch, switch) 
    
    model = abm.ABM(Voter{nissues}, space;properties = properties)

    add_agents!(model, nagents, nissues, voter_pos_initializor)

    set_partyposs(model, party_pos_hardwired)

    update_pposs_tracker!(model)

    update_voters_partyids!(model)

    initialize_vpids_tracker(model)

    initialize_incumbent_streak_counter!(model)

    initialize_ballot_trackers(model)
    
    initialize_median_pos_tracker(model)

    return(model)
end


initialize_model(x::ModelParams) = initialize_model(;ntfromstruct(x)...)


# * Stepping



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

candidates_iteration_setup!(m) = candidates_iteration_setup!(m, m.properties[:switch])

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

function get_mean_among_supporters(supporters::Vector,
     model,
      weights)
        [mean([model[i].pos[issue]
              for i in supporters], weights)
         for issue in 1:model.properties[:nissues]]
end



#= 
function get_new_parties_poss(model, new_supporters, old_supporters)
    ω = model.properties[:ω]

    mean_previous_supporters = dictmap(v->get_mean_among_supporters(v,model),
                                       old_supporters )

    mean_new_supporters = kvdictmap((k,v)-> isempty(v) ? model[k].pos : get_mean_among_supporters(v,model),
                                  new_supporters)

    kvdictmap((k,v)-> (ω .* v .+ ((1-ω) .* mean_new_supporters[k])) |> Tuple, mean_previous_supporters )

end
=#
function set_new_parties_poss!(model,newpposs)
    for (k,v) in newpposs
        model[k].pos = v
        model.properties[:parties_candidateid_ppos_δ][k][:partyposition] = v
    end
end 

function get_new_parties_poss(model)
    supporters =  get_parties_supporters(model)
    loyalties(listofagents) = map(a -> loyalty(a,model), listofagents )
    weightsdict = dictmap(l -> Weights(loyalties(l)), supporters)
    newppossdict = kvdictmap((k,v)->Tuple(get_mean_among_supporters(v,model, weightsdict[k])),
    supporters )
    
    return(newppossdict)
end

function set_agent_new_κ!(agentid,model)
    myLast_PartyVote = model.properties[:voterBallotTracker][agentid][end]
    proportion_IvotedForThisParty = proportionmap(model.properties[:voterBallotTracker][agentid])[myLast_PartyVote]
    baseline_κ = model.properties[:κ]
    model[agentid].κ = proportion_IvotedForThisParty * baseline_κ
end

# FIXME: this is a loop per step for stuff I'm not even using lol
function set_agents_new_κ!(model)
    kappa_switch = model.properties[:kappa_switch]
    for i in abm.allids(model)
        if kappa_switch == :off
            break
        else
            set_agent_new_κ!(i, model)
        end
    end
end

# this function must be run AFTER
# the candidates_iteration_setup!
function pre_stepping_loop!(model)
       for i in abm.allids(model)
            candidates = map(x->x[:partycandidate],
            values(model.properties[:parties_candidateid_ppos_δ]))
            myvote = get_closest_fromList(i,candidates,model)
            push!(model.properties[:voterBallotTracker][i], model[myvote].myPartyId)
        end
end

function first_step_loop!(model)
    candidates_iteration_setup!(model, :initial)
    pre_stepping_loop!(model)
end

function second_third_steps_loop!(model)
    candidates_iteration_setup!(model, :second)
    pre_stepping_loop!(model)
end

function turn_currrent_incumbent_into_old!(model)
        model.properties[:incumbent_streak_counter].old_incumbentholder = model.properties[:incumbent_party]
end

function get_new_winner_party(model)
    new_winner_party = model[getmostvoted(model, :iteration)].myPartyId
    return(new_winner_party)
end

function update_incumbent!(model)
    new_winner = get_new_winner_party(model)
    model.properties[:incumbent_party] = new_winner
end

function update_voters_ballots_history!(model)
        for i in abm.allids(model)
            party_i_votedfor = get_whichCandidatePartyAgentVotesfor(i,model)[2]
            push!(model.properties[:voterBallotTracker][i], party_i_votedfor )
        end
end


function update_within_party_shares!(model)
    model.properties[:withinpartyshares] = get_withinpartyshares(model)
end


function update_cross_voting_tracker!(model)
    push!(model.properties[:cross_voting], 0)
    for i in abm.allids(model)
        add_crossvoting_tocounter!(i, model)
    end
end

function  update_keep_probs_tracker!(model)
    nagents = abm.nagents(model)
    if length(model.properties[:keep_probs]) != nagents
        model.properties[:keep_probs] = Vector(undef,nagents)
    end
                
    for i in abm.allids(model)
        model.properties[:keep_probs][i] =  get_keep_party_id_prob(i,model)
    end
end

function update_partyid_tracker!(i, model)
    model.properties[:voters_partyids][i] = model[i].myPartyId
end

function update_partyids!(model)
    push!(model.properties[:party_switches], 0)
    # this function can only be run AFTER update_keep_probs!
    for i in abm.allids(model)
        keep_prob = model.properties[:keep_probs][i]
        update_partyid!(i,model,keep_prob)
        update_partyid_tracker!(i,model)
    end
end 


function update_new_parties_poss!(model)
    newposs = get_new_parties_poss(model)
    set_new_parties_poss!(model,newposs)
end

function assume_initial_partyid!(model)

    for i in abm.allids(model) assume_initial_partyid!(i, model)
    end
end

function model_actual_loop!(model)
        assume_initial_partyid!(model)
        candidates_iteration_setup!(model)
        set_agents_new_κ!(model)
        turn_currrent_incumbent_into_old!(model)
        update_incumbent!(model)
        update_streakCounter!(model)
        update_voters_ballots_history!(model)
        update_within_party_shares!(model)
        update_cross_voting_tracker!(model)
        update_keep_probs_tracker!(model)
        update_partyids!(model)
        update_new_parties_poss!(model)
end

function model_step!(model)
    model.properties[:is_at_step] += 1
    if model.properties[:is_at_step] == 1
        first_step_loop!(model)
    elseif (model.properties[:is_at_step] == 2) || (model.properties[:is_at_step] == 3)
        second_third_steps_loop!(model)
    else
        model_actual_loop!(model)
    end
end

foursteps!(m) =  for _ in 1:4 model_step!(m) end
