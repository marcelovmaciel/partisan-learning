import Pkg
Pkg.activate(".")

import Dictionaries as DICT
import Agents as abm
import Distributions as distri
import Distances as dist

using AlgebraOfGraphics, GLMakie
using Deldir
using StaticArrays
using InteractiveDynamics

# * The initial condition

# ** Model Attributes
v = 3000 # stub
p = 10 # stub
n = 2
k = 0:8
s = 0:2
c = 10;

properties = Dict(:r => 3, :c => 10, :partyids => Int64[])


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

function model(;
               n = n,
               k = k,
               s = s,
               v = v,
               p = p,
               properties = properties)
    numagents = v + p
    space = abm.ContinuousSpace((last(k), last(k)))
    model = abm.ABM(Union{AnchoredParty,Voter},space, properties = properties)

    for i in 1:v
        agent = Voter(i,n,k,s)
        abm.add_agent!(agent,model)
    end


   for i in 1:v
        model[i].id_pos = MVector{n}(model[i].pos...)

     #   println(model[j])
    end

    dummy_share = 0.

    for j in v+1:((v+p))
        agent = AnchoredParty(j,n,k,dummy_share)
        abm.add_agent!(agent,model)
    #    println(agent)
        end

    for j in v+1:((v+p))
        model[j].id_pos = MVector{n}(model[j].pos...)

     #   println(model[j])
    end

    partyids = collect(v+1:(v+p))
    model.properties[:partyids] = partyids
#    println(model[v+1])
    set_shares!(model)

    # this is potentially wrong!!! The radius_share differs from the whole share!!!
    for j in v+1:((v+p))
        model[j].potential_pos = (position = model[j].id_pos,  radius_share = model[j].share)
     #   println(model[j])
    end

    return model

end

# *  The stepping
# ** Salience signal
function get_voters_withinradius(p,voters)::Vector{Voter}
    voters[(a->dist.euclidean(a.id_pos,p.id_pos)).(voters) .<  m.properties[:r]]
end

function getparties_spheres(parties,voters)
    map(p->get_voters_withinradius(p, voters), parties)
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

function positive_signalfn(x::Voter, positive_signal_index)
    weight = copy(x.weight)
    if weight[positive_signal_index] == last(s)
       weight[positive_signal_index] += 0
    else
        weight[positive_signal_index] += 1
    end
        return(weight)
end

function negative_signalfn(x::Voter, nsi)
    weight = copy(x.weight)
    if weight[nsi] == first(s)
        weight[nsi] += 0
    else
        weight[nsi] -= 1
    end
    return(weight)
end

function shares_psphere(p, prties_spheres, parties)
    share = vote_shares(m.properties[:partyids], map(a->closest_party(a,parties), prties_spheres[p.id]))
end

function closest_party_newWeight(agent, new_weight, parties = parties)
    oldweight = copy(agent.weight)
    agent.weight = new_weight
    cp = closest_party(agent, parties)
    agent.weight = oldweight
    return(cp)
    end

function testsignals(p, prties_spheres; n=n, parties = parties)

    signal = whichsignal(n)

    test_positive = vote_shares(m.properties[:partyids],
                                map(x-> closest_party_newWeight(x, positive_signalfn(x, signal.positive), parties), prties_spheres[p.id]))

    test_negative = vote_shares(m.properties[:partyids],
                                map(x-> closest_party_newWeight(x, negative_signalfn(x, signal.negative), parties), prties_spheres[p.id]))
    return(positive = test_positive, negative = test_negative, signal = signal)
end

function send_signal(p, prties_spheres = prties_spheres, parties = parties)


    status_quo_share = shares_psphere(p, prties_spheres, parties)[p.id]

    counterfactuals = testsignals(p, prties_spheres, parties = parties);
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

function sampleposition_withinradius(p,r=m.properties[:r], k = k)
    x,y = p.pos

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
                          voters = voters)

    status_quo_share = shares_psphere(p, prties_spheres,parties)[p.id]
#    println(status_quo_share))

    old_position = copy(p.id_pos)

    new_position = sampleposition_withinradius(p)

    p.id_pos = new_position

    newpartieslist = DICT.Dictionary(m.properties[:partyids],
    map(i->m[i],m.properties[:partyids]))

    newprties_spheres = getparties_spheres(newpartieslist, voters)

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
function update_weight!(voter, party, prties_spheres = prties_spheres )

    similarity = 1-dist.cosine_dist(voter.id_pos, party.id_pos)
    mc = distri.Uniform(0,1)
    if rand(mc) < similarity && party.signal.shouldsend.positive
        voter.weight = positive_signalfn(voter, party.signal.whichsignal.positive)
        end

    if rand(mc) < similarity && party.signal.shouldsend.negative
        voter.weight = negative_signalfn(voter,party.signal.whichsignal.negative)
    end
end

# ** Opinion Dynamics

function axelrod_update!(voter, model; n = n)

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
m = model()

function signalstep!(m, parties, voters, prties_spheres)

    for p in parties
        p.signal = send_signal(p, prties_spheres, parties)
    end

    for v in voters
        for p in parties
            update_weight!(v,p,prties_spheres)
        end
    end
end


function opinionstep!(m; voterid = 1:v)
    for i in voterid
        axelrod_update!(m[i],m)
    end
    for i in voterid
        setfield!(m[i], :pos , Tuple(m[i].id_pos))
    end
end

function pupdate_potential_pos!(p, parties, prties_spheres,voters)

    position, radius_share = test_newposition(p, prties_spheres,  parties, voters)

    if radius_share > p.potential_pos.radius_share
        p.potential_pos = (position=position, radius_share = radius_share)
    end

end

function allupdate_potential_pos!(m, parties,voters,
                                  prties_spheres)
    partyids = m.properties[:partyids]
    for p in partyids
        pupdate_potential_pos!(m[p], parties, prties_spheres, voters)
    end
end

function campaign_cycle!(m; c = c)

    parties = getparties(m)
    voters = let vs = Vector{Voter}(undef,v)
        for i in 1:v
            vs[i] = (@inbounds m[i])
        end
        vs
    end

    prties_spheres = getparties_spheres(parties, voters);

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

#m = model()

#abm.run!(m, abm.dummystep, model_step!, 10)

#m[3001] |> typeof |> fieldnames


# * Visualization



agent_colors(a) = typeof(a) == Voter{2} ? "#2b2b33" :  "#bf2642"

agent_size(a) = typeof(a) == Voter{2} ? 5 :  20
m = model()

# fig, abmstepper = abm_plot(m; ac = agent_colors, as = agent_size)
# fig




# abm.run!(m, abm.dummystep, model_step!, 10)

# fig2, abmstepper2 = abm_plot(m; ac = agent_colors)
# fig2



m = model()
abm_play(
    m,
    abm.dummystep,
    model_step!,
    ac = agent_colors,
    as = agent_size
)




parties = (collect(values(getparties(m))))

_, vor, _ = deldir(
    map(x-> x.id_pos[1],parties),
       map(x-> x.id_pos[2],parties),
       [0.,8.1,0.0,8.1]);

Vx, Vy = edges(vor);

data1 = map(x-> x.id_pos[1],abm.allagents(m));
data2 = map(x-> x.id_pos[2],abm.allagents(m));

f = Figure(resolution= (400,200))

scatter(f[1,1], data1, data2, markersize = 2)
scatter!(f[1,1], map(x-> x.id_pos[1],parties ),
         map(x-> x.id_pos[2],parties ), color = :red)
lines!(f[1,1],Vx,Vy)

Makie.save("firstplot.png",f)


# ** See weights update
resolution = (400, 400)
fig = Figure(; resolution)
ax = Axis(fig[1, 1], title="Some plot")
ax2 = Axis(fig[1, 2], title="Some plot")

foo = map(x-> x.weight,filter(x->typeof(x)== Voter{2}, collect(abm.allagents(m))))

hist!(fig[1,1], map(x->x[1], foo))
hist!(fig[1,2], map(x->x[2], foo))
fig

signalstepping(m)




is, ss = collect(zip([(k,v) for (k,v) in pairs(shares)]...));
is = [i for i in is];
ss = [s for s in ss];
df = (votes = [Int64(i) for i in filter(x->!isnothing(x),votes)], is=is)

# for ae in ag
#     Axis(ae).xticklabelrotation[] = π/2
# end
