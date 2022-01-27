import Pkg

Pkg.activate("../")

using Agents
using Distributions
using GLMakie
using InteractiveDynamics

mutable struct DummyAgent{n} <: AbstractAgent
    id::Int
    pos::NTuple{n,Float64}
    am_I_plus::Bool
    var_to_plot::Float64
end


function initialize_model(nagents=500, n=2)

    space = ContinuousSpace(ntuple(x -> float(10),n), periodic = false)
    model = ABM(DummyAgent{n}, space)
    for i in 1:nagents
        pos = Tuple(rand(Uniform(0,10), n))
        add_agent_pos!(DummyAgent{n}(i, pos, false, rand(Uniform(0,10))), model)
    end

    for i in rand(allids(model), 250)
        model[i].am_I_plus = true
    end
    return(model)
end

function var_increase(a)
    a.am_I_plus ? a.var_to_plot += 1 : a.var_to_plot -= 1
end


function model_step!(model)
    for i in allids(model)
        var_increase(model[i])
    end
end

m = initialize_model()


fig,adf,mdf = abm_data_exploration(m,
                                     dummystep,
                                     model_step!)

var_wanna_observe = [Observable([m[i].var_to_plot * 10]) for i in  allids(m)]
nsteps = Observable([0.])
function newstep(m, var_wanna_observe = var_wanna_observe, nsteps = nsteps)
      model_step!(m)
      var_values = [m[i].var_to_plot * 10 for i in  allids(m)]
      for (i,v) in enumerate(var_values)
          var_wanna_observe[i][] = [v]
      end
      nsteps[] = push!(nsteps.val,nsteps.val[end]+ 1.)
  end

fig,adf,mdf = abm_data_exploration(m, dummystep, newstep)



scatter(fig[1,2], nsteps, var_wanna_observe[1])
lines!(fig[1,2],nsteps, var_wanna_observe[1])
for i in var_wanna_observe[2:end]
    scatter!(fig[1,2], nsteps, i)
    lines!(fig[1,2], nsteps, i)
end
