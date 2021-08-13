import Pkg
Pkg.activate("../")

using Revise

using PartisanLearning
# * Visualization

foo = ModelParams()

foo |> typeof

Parameters.@unpack foo

param.v,
param.p,
param.n,
param.k,
param.s,
param.c,
param.r,
param.voterids,
param.partyids,

v, p, n, k, s, c, r, voterids, partyids = Parameters.@unpack foo

Dict((x->(fn=>getfield(x, fn) for fn ∈ fieldnames(typeof(x))))(foo))

typeof(foo) |>fieldnames

typedict(foo)


agent_colors(a) = typeof(a) == Voter{2} ? "#2b2b33" :  "#bf2642"

agent_size(a) = typeof(a) == Voter{2} ? 5 :  20


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
