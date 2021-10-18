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
# * Issue Salience Model
fooparams = is.IssueSalience.ModelParams()

m =  is.model(fooparams)

is.model_step!(m)

m.properties

agent_colors(a) = typeof(a) == is.Voter{2} ? "#2b2b33" :  "#bf2642"

agent_size(a) = typeof(a) == is.Voter{2} ? 5 :  20

m = model()
pl.abm_play(
    m,
    is.abm.dummystep,
    is.model_step!,
    ac = agent_colors,
    as = agent_size
)


# * Unorganized code


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

foo = map(x-> x.weight,filter(x->typeof(x)== Voter{2},
collect(abm.allagents(m))))

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
# e
end



# * PartyidModel



m.properties[:partyids]

ids = collect(allids(m))

parties = pla.sample_parties_pos(2,m)

candidate = pid.StatsBase.sample(collect(nearby_ids(m[parties[1]].pos,m,0.1, exact = true)))

a = Tuple[]

push!(a, (1,1))

a

collect((pid.StatsBase.sample(collect(nearby_ids(p,m,0.1, exact = true))) for p in parties))

pid.dist.euclidean(m[parties[1]].pos, m[candidate].pos)




agent_colors(a) = a.id == m.properties[:incumbent]  ? :yellow : (a.amIaCandidate  ?  "#bf2642"  : "#2b2b33")
agent_size(a) = a.id == m.properties[:incumbent]  ? 20 : (a.amIaCandidate ? 15 : 5)

# fig |> typeof

params = Dict(:κ => 0.0:0.1:1,
              :ρ => 0.0:0.1:1,
              :ϕ => 0.0:0.1:1)


# So, bigger ρ means agents will converge faster to the mean
# partyid
# higher phi has the same effect
#higher kappa means agents will tend to vote more for
#their partyid


adata = [(a->(pid.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
         (a->(pid.get_distance_IvsParty(a,m)), pid.StatsBase.mean)]

alabels = ["Voted Against PartyId", "dist(closest,party's)"]

fig,adf,mdf = abm_data_exploration(m,
                                   pid.abm.dummystep,
                                   pid.model_step!,
                                   params;
                                   adata, pid.mdata,
                                   alabels,
                                   ac = agent_colors,
                                   as = agent_size, spu = 1)

# ** Trying to add the partyid stuff
#
 fig |> typeof |> fieldnames

ncandidates = 3
nissues = 2
m = pid.initialize_model(1000,nissues, ncandidates)

agent_colors(a) = a.id == m.properties[:incumbent]  ? :yellow : (a.amIaCandidate  ?  "#bf2642"  : "#2b2b33")
agent_size(a) = a.id == m.properties[:incumbent]  ? 20 : (a.amIaCandidate ? 15 : 5)

fig,stepper = pl.abm_play(
    m,
    pid.abm.dummystep,
    pid.model_step!,
    ac = agent_colors,
    as = agent_size,
)

# fig |> typeof |> fieldnames
# fig.scene |> typeof |> fieldnames

#pids = (model) -> ([Point{2}([m[x].myPartyId[1],m[x].myPartyId[2]]) for x in  pid.abm.allids(model)])
#collect(Point(x) for x in values(m.properties[:partyids]))
#foo = Observable(collect(values(m.properties[:partyids])))

#scatter(fig[1,2], foo)

# fig


#foo[] = collect(values(m.properties[:partyids]))


#fig
# fig2
foo = Observable(collect(values(m.properties[:partyids])))


scatter(fig[1,2], foo)

function newstep(m, foo = foo)
    pid.model_step!(m)
    foo[] = collect(values(m.properties[:partyids]))
end


function abm_play!(fig, abmstepper, model, agent_step!, model_step!; spu,
                   newvartoplot, obsstepper)
    # preinitialize a bunch of stuff
    model0 = deepcopy(model)
    modelobs = Observable(model) # only useful for resetting
    speed, slep, step, run, reset, = InteractiveDynamics.abm_controls_play!(fig, model, spu, false)

    #lift(x->newplotter(fig[newviewerpos...], x), collect(values(myobs.val))) # Add it to the figure


# fig
    # Clicking the step button
    on(step) do clicks
        n = speed[]
        Agents.step!(abmstepper, model, agent_step!, model_step!, n)
    end

    # Clicking the run button
    isrunning = Observable(false)
    on(run) do clicks; isrunning[] = !isrunning[]; end
    on(run) do clicks
        @async while isrunning[]
            n = speed[]
            model = modelobs[] # this is useful only for the reset button
            Agents.step!(abmstepper, model, agent_step!, model_step!, n)


            slep[] == 0 ? yield() : sleep(slep[])
            isopen(fig.scene) || break # crucial, ensures computations stop if closed window.
        end
    end
    # Clicking the reset button
    on(reset) do clicks
        modelobs[] = deepcopy(model0)
        Agents.step!(abmstepper, modelobs[], agent_step!, model_step!, 0)
    end
    return nothing
end

function obsstep(m)
    collect(values(m.properties[:partyids]))
end

abm_play!(fig,stepper, m, pid.abm.dummystep,
          newstep; newvartoplot = foo,
          obsstepper = obsstep, spu = 0.1)
fig
# function static_preplot!(ax, model)
#     obj = CairoMakie.scatter!([50 50]; color = :red) # Show position of teacher
#     CairoMakie.hidedecorations!(ax) # hide tick labels etc.
#     CairoMakie.translate!(obj, 0, 0, 5) # be sure that the teacher will be above students
# end
# idxes = collect(pid.abm.allids(m))
# xs = Observable([m[x].myPartyId[1] for x in  idxes])
# ys = Observable([m[x].myPartyId[2] for x in  idxes])


fig


# oxs = on(xs) do val
#     val
# end

# oys =  on(ys) do val
#     val
# end

ax, scat = scatter(fig[1,3],xs,ys)




fig.layout

fig2,stepper = pl.abm_play(
    m,
    pid.abm.dummystep,
    pid.model_step!,
    ac = agent_colors,
    as = agent_size,
)


stepper.pos

scatter(fig2[1,3],xs,ys)


fig |> typeof |> fieldnames

fig.content[end-1]

on





# * PartyLabelModel
ncandidates = 10
nissues = 5

# BUG: run the δ thing might be leading to a bug!
m = pla.initialize_model(1000,nissues,
                         ncandidates)
m.properties[:incumbent]

pla.candidates_iteration_setup!(m)

m.properties[:partiesposs]


proportion_peers_voteLikeMe2(1,m)

withinpartyshares[m.properties[:voters_partyids][1]]



m.properties[:voterBallotTracker]

filter((t->t[2]==810),collect(m.properties[:voters_partyids]))

collect(m.properties[:voters_partyids])

m.properties[:voters_partyids]




sample(collect(nearby_ids(m[208],
                       m,
                       0.3,
                       exact = true)))

m.properties[:partiesposs][208]

for (k,v) in m.properties[:partiesposs]
println(k, " ", v[:partycandidate])
end


m.properties[:voterBallotTracker]

m.properties



for i in allids(m)
    println(pla.get_closest_candidate(i,m))
end


for k in keys(m.properties[:partiesposs])
    println(k, " ", m.properties[:partiesposs][k][:partycandidate] )
end


# * Test this
using Agents
using InteractiveDynamics
using GLMakie
using Statistics

model, agent_step!, model_step! = Models.schelling()

params = Dict(:min_to_be_happy => 1:1:5)

fig, adf, mdf = abm_data_exploration(model, agent_step!, model_step!, params;
                                     adata=[(:mood, mean)])


count_unhappy(model) = count(a.mood == false for a in allagents(model))
fig, adf, mdf = abm_data_exploration(model, agent_step!, model_step!, params;
    adata=[(:mood, mean)], mdata=[count_unhappy])
