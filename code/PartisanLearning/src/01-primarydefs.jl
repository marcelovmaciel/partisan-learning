# * Issue Salience model

module IssueSalience

import Dictionaries as DICT
import Agents as abm
import Distributions as distri
import Distances as dist
import Base.@kwdef
using StaticArrays

# * The initial condition
# ** Model Attributes


@kwdef struct ModelParams
    v = 3000
    p = 10
    n = 2
    k = 0:8
    s = 0:2
    c = 10
    r = 3
    voterids = 1:v
    partyids = v+1:((v+p))
end


# ** Voters Attributes
mutable struct Voter{n} <: abm.AbstractAgent
    id::Int
    pos::NTuple{n,Float64}
    id_pos::MVector{n,Float64}
    weight::MVector{n,Int}
end

function voter_issue_generation(n,k)
    # anchor = rand(distri.DiscreteUniform(first(k), last(k)))
    anchor = rand(distri.Uniform(first(k), last(k)))
    function anchored_positions(anchor)
        # anchorp1 = (anchor == last(k) ? anchor : anchor + 1 )
        # anchorm1 = (anchor == first(k) ? anchor : anchor - 1 )
        anchorp1 = (anchor + 1 > last(k) ? last(k) : anchor + 1 )
        anchorm1 = (anchor - 1 < first(k) ? first(k) : anchor - 1 )
        return(anchorm1, anchorp1)
        end
    MVector{n}(rand(distri.Uniform(anchored_positions(anchor)...), n,1)) #Tuple(rand(anchored_positions(anchor), n,1))
end

function weight_assigner(n,s)
    MVector{n}(rand(collect(s), n,1))
end

function Voter(id,n,k,s)
    ideology= voter_issue_generation(n,k)
    Voter{n}(id,
            Tuple(ideology),
             ideology,
             weight_assigner(n,s))

end

# ** Parties attributes

mutable struct AnchoredParty{n} <: abm.AbstractAgent
    id::Int
    pos::NTuple{n,Float64}
    id_pos::MVector{n,Float64}
    share::Float64
    signal::NamedTuple{(:shouldsend, :whichsignal),
                       Tuple{NamedTuple{(:positive, :negative),
                                        Tuple{Bool, Bool}},
                             NamedTuple{(:positive, :negative),
                                        Tuple{Int64, Int64}}}}
    potential_pos::NamedTuple{(:position, :radius_share), Tuple{MVector{n, Float64}, Float64}}
end

mutable struct UniformParty{n} <: abm.AbstractAgent
    id::Int
    pos::MVector{n,Float64}
    id_pos::MVector{n,Float64}
    share::Float64
end


function AnchoredParty(id,n,k,share)
    generate_party_ideology(n,k) = voter_issue_generation(n,k)
    ideology = generate_party_ideology(n,k)
    dummysignal = (shouldsend = (positive = true, negative = true),
                   whichsignal = (positive = 1, negative = 1))
    AnchoredParty{n}(id, Tuple(ideology), ideology, 0., dummysignal, (position = ideology, radius_share = 0.))

 end


function UniformParty(id,n,k,share)
    generate_party_ideology(n,k) =  MVector{n}(rand(distri.Uniform(first(k), last(k)),n,1))
    ideology = generate_party_ideology(n,k)
    UniformParty{n}(id, Tuple(ideology), ideology, 0.)
end


# ** Secondary attribute fillers

utility(agent,party) = -dist.weuclidean(agent.id_pos, party.id_pos, agent.weight)

function getparties(model)
    parties = DICT.Dictionary(model.properties[:partyids],
                              map(x->model[x], model.properties[:partyids]))

end


function closest_party(agent::Voter, parties)
    argmax(Dict(map(p-> p.id=>utility(agent, p), parties)))
end


closest_party(agent::AnchoredParty, parties) = nothing

function vote_shares(partyid, votes)
    d = Dict(zip(partyid, zeros(length(partyid))))
    for val in partyid
        for vote in votes
            if isa(vote, Number) && isnan(vote)
            continue
        end
            if vote == val
                d[val] = get(d, val, 0) + 1
            else
                d[val] += 0
            end
        end
    end

    s = DICT.sum(DICT.Dictionary(d))

    lastd = map(x->x/s, DICT.Dictionary(d))

    return lastd
end

function set_shares!(model)
    parties = DICT.Dictionary(model.properties[:partyids],
                              map(x->model[x], model.properties[:partyids]))
    votes = map(a->closest_party(a,parties), abm.allagents(model));
    shares = vote_shares(model.properties[:partyids], votes)
    for agent_index in model.properties[:partyids]
        model[agent_index].share = shares[agent_index]
    #    println(m[agent_index].pos,m[agent_index].id_pos)
    end
end

# ** Model initialization

function model(param::ModelParams)

    v, p, n, k, s, c, r, voterids, partyids = param.v, param.p, param.n, param.k, param.s, param.c, param.r, param.voterids, param.partyids
    function typedict(params)
        Dict((x->(fn=>getfield(x, fn) for fn ∈
            fieldnames(typeof(x))))(param))
    end

    properties = typedict(param)

    space = abm.ContinuousSpace((last(k), last(k)))
    model = abm.ABM(Union{AnchoredParty,Voter},space, properties = properties)

    for i in voterids
        agent = Voter(i,n,k,s)
        abm.add_agent!(agent,model)
    end

    for i in voterids
        model[i].id_pos = MVector{n}(model[i].pos...)
    end

    dummy_share = 0.

    for j in partyids
        agent = AnchoredParty(j,n,k,dummy_share)
        abm.add_agent!(agent,model)
        end
    for j in partyids
        model[j].id_pos = MVector{n}(model[j].pos...)
    end
    set_shares!(model)

    # this is potentially wrong!!! The radius_share differs from the whole share!!!
    for j in partyids
        model[j].potential_pos = (position = model[j].id_pos,  radius_share = model[j].share)
    end
    return model
end


# *  The stepping
# ** Salience signal
function get_voters_withinradius(p,voters, m)::Vector{Voter}
    voters[(a->dist.euclidean(a.id_pos,p.id_pos)).(voters) .<  m.properties[:r]]
end

function getparties_spheres(parties,voters,m)
    map(p->get_voters_withinradius(p, voters,m), parties)
end

function whichsignal(nissues)
    issues = collect(1:nissues)
    pos1 = rand(issues)
    pos2 = rand(issues)
    if pos1 == pos2
        deleteat!(issues, findall(x->x==pos1,issues)) # deletes all elements that have a certain value from an array
        pos2 = rand(issues)
        end
    signals = (positive = pos1 , negative = pos2 )
    end

function positive_signalfn(x::Voter, positive_signal_index;m=m)
    weight = copy(x.weight)
    if weight[positive_signal_index] == last(m.properties[:s])
       weight[positive_signal_index] += 0
    else
        weight[positive_signal_index] += 1
    end
        return(weight)
end

function negative_signalfn(x::Voter, nsi;m=m)
    weight = copy(x.weight)
    if weight[nsi] == first(m.properties[:s])
        weight[nsi] += 0
    else
        weight[nsi] -= 1
    end
    return(weight)
end

function shares_psphere(p, prties_spheres, parties,m)
    share = vote_shares(m.properties[:partyids], map(a->closest_party(a,parties), prties_spheres[p.id]))
end

function closest_party_newWeight(agent, new_weight, parties = parties)
    oldweight = copy(agent.weight)
    agent.weight = new_weight
    cp = closest_party(agent, parties)
    agent.weight = oldweight
    return(cp)
    end

function testsignals(p, prties_spheres; parties = parties,m=m)

    signal = whichsignal(m.properties[:n])

    test_positive = vote_shares(m.properties[:partyids],
                                map(x-> closest_party_newWeight(x, positive_signalfn(x, signal.positive,m=m), parties), prties_spheres[p.id]))

    test_negative = vote_shares(m.properties[:partyids],
                                map(x-> closest_party_newWeight(x, negative_signalfn(x, signal.negative,m=m), parties), prties_spheres[p.id]))
    return(positive = test_positive, negative = test_negative, signal = signal)
end

function send_signal(p, prties_spheres = prties_spheres, parties = parties, m=m)


    status_quo_share = shares_psphere(p, prties_spheres, parties,m)[p.id]

    counterfactuals = testsignals(p, prties_spheres, parties = parties,m=m);
        positive_counterfactual_share = counterfactuals.positive[p.id]
    negative_counterfactual_share = counterfactuals.negative[p.id]
    send_positive = positive_counterfactual_share > status_quo_share
    send_negative = negative_counterfactual_share > status_quo_share

    shouldsend = (positive= send_positive,
                  negative = send_negative)
    whichsignal = counterfactuals.signal
    return(shouldsend = shouldsend, whichsignal=whichsignal)
end


# ** Explorer

function sampleposition_withinradius(p, m)
    x,y = p.pos
    r = m.properties[:r]
    k = m.properties[:k]
    function cropvalue(tentativevalue)
        if tentativevalue > last(k)
            tentativevalue = Float64(last(k))
        elseif tentativevalue < first(k)
            tentativevalue = Float64(first(k))
        else
            tentativevalue = tentativevalue
        end
        return(tentativevalue)
    end

    newx = cropvalue(rand(distri.Uniform(x-r, x+r)))
    newy = cropvalue(rand(distri.Uniform(y-r,y+r)))

    return(newx,newy)
end


function test_newposition(p,
                          prties_spheres = prties_spheres,
                          parties= parties,
                          voters = voters,
                          m=m)

    status_quo_share = shares_psphere(p, prties_spheres,parties,m)[p.id]
#    println(status_quo_share))

    old_position = copy(p.id_pos)

    new_position = sampleposition_withinradius(p,m)

    p.id_pos = new_position

    newpartieslist = DICT.Dictionary(m.properties[:partyids],
    map(i->m[i],m.properties[:partyids]))

    newprties_spheres = getparties_spheres(newpartieslist, voters,m)

    new_share = vote_shares(
        m.properties[:partyids],
        map(a->closest_party(a, newpartieslist), newprties_spheres[p.id]))[p.id]

    p.id_pos = old_position

    share = 0.


    if new_share > status_quo_share
        ret = new_position
        share += new_share
    else
        ret = old_position
        share += status_quo_share
        end
    return(position = ret, radius_share = share)
end


# ** Voter weight update
function update_weight!(voter, party, prties_spheres = prties_spheres,m=m )

    similarity = 1-dist.cosine_dist(voter.id_pos, party.id_pos)
    mc = distri.Uniform(0,1)
    if rand(mc) < similarity && party.signal.shouldsend.positive
        voter.weight = positive_signalfn(voter, party.signal.whichsignal.positive,m=m)
        end

    if rand(mc) < similarity && party.signal.shouldsend.negative
        voter.weight = negative_signalfn(voter,party.signal.whichsignal.negative,m=m)
    end
end

# ** Opinion Dynamics

function axelrod_update!(voter, model)
    n = model.properties[:n]
    neighbor = abm.random_agent(model, x-> typeof(x) == Voter{n})

    voterpos,neighborpos = MVector{n}(voter.pos),MVector{n}(neighbor.pos)


    similarity = 1-dist.cosine_dist(voterpos, neighborpos)
    mostsalient_issue = argmax(voter.weight)

    mc = distri.Uniform(0,1)

    function newissuepos()
        ((voter.id_pos[mostsalient_issue]*voter.weight[mostsalient_issue] +
            neighbor.id_pos[mostsalient_issue]*neighbor.weight[mostsalient_issue])/
            (neighbor.weight[mostsalient_issue] +
            voter.weight[mostsalient_issue]))
    end
    if rand(mc) < similarity
        voter.id_pos[mostsalient_issue] = newissuepos()
    end
end


# ** The actual stepping

function signalstep!(m, parties, voters, prties_spheres)

    for p in parties
        p.signal = send_signal(p, prties_spheres, parties,m)
    end

    for v in voters
        for p in parties
            update_weight!(v,p,prties_spheres,m)
        end
    end
end


function opinionstep!(m)

    for i in m.properties[:voterids]
        axelrod_update!(m[i],m)
    end

    for i in m.properties[:voterids]
        setfield!(m[i], :pos , Tuple(m[i].id_pos))
    end
end

function pupdate_potential_pos!(p, parties, prties_spheres,voters,m)

    position, radius_share = test_newposition(p, prties_spheres,  parties, voters,m)

    if radius_share > p.potential_pos.radius_share
        p.potential_pos = (position=position, radius_share = radius_share)
    end

end

function allupdate_potential_pos!(m, parties,voters,
                                  prties_spheres)
    partyids = m.properties[:partyids]
    for p in partyids
        pupdate_potential_pos!(m[p], parties, prties_spheres, voters,m)
    end
end

function campaign_cycle!(m)
    c = m.properties[:c]

    parties = getparties(m)
    voters = let vs = Vector{Voter}(undef,m.properties[:v])
    for i in m.properties[:voterids]
            vs[i] = (@inbounds m[i])
        end
        vs
    end

    prties_spheres = getparties_spheres(parties, voters,m);

    for _ in c
        signalstep!(m, parties, voters, prties_spheres)
        opinionstep!(m)
        allupdate_potential_pos!(m,parties,voters,prties_spheres)
    end

end


function electoral_iteration!(m)
    partyids = m.properties[:partyids]
    for p in partyids
        m[p].pos = Tuple(m[p].potential_pos.position)
        m[p].id_pos = m[p].potential_pos.position
    end
    set_shares!(m)
end

function model_step!(m)
    campaign_cycle!(m)
    electoral_iteration!(m)
 end


# My guess is that parties are seeing the OLD voters positions, rather than
# their actual position!
#
end


# * Party id model

module PartyId


end
