using Test
import PartisanLearning as pl
const is = pl.IssueSalience
const pid = pl.PartyId



@testset "Is the constructor well formed?" begin
    v1 = let amIaCandidate = false
        nissues = 1
        pos = Tuple(1.)
        myPartyId = pos
        return(pid.Voter{nissues}(1,pos,amIaCandidate, myPartyId))
    end
    @test v1 == pid.Voter(1,1,Tuple(1.))
end


function getstructvalues(S)
    map(s -> getfield(S,s), fieldnames(typeof(S)))
end

@testset "Does the package change my agent's position?" begin
    v2 = let space = pid.abm.ContinuousSpace((2.,))
        nissues = 1
        model = pid.abm.ABM(pid.Voter{nissues}, space)
        pid.abm.add_agent_pos!(pid.Voter(1, nissues, Tuple(1.,)), model)
        pid.abm.allagents(model) |> first
    end
    v3 = pid.Voter(1,1,Tuple(1.,))

    @test getstructvalues(v2) == getstructvalues(v3)
    v4 = let space = pid.abm.ContinuousSpace((2.,))
        nissues = 1
        model = pid.abm.ABM(pid.Voter{nissues}, space)
        pid.abm.add_agent!(Tuple(1.), model,false, (1.,))
        pid.abm.allagents(model) |> first
    end
    @test  getstructvalues(v4) == getstructvalues(v2)
end


#= So, I do not generally want to use the add_agent! procedure, because i tend
to control my agents' constructors =#

function one_agent_model()
    space = pid.abm.ContinuousSpace((2.,))
    nissues = 1
    model = pid.abm.ABM(pid.Voter{nissues}, space)
    pid.abm.add_agent!(Tuple(1.), model,false, (1.,))
    model
end

function simple_agents_model(nagents)
    space = pid.abm.ContinuousSpace((2.,))
    nissues = 1
    model = pid.abm.ABM(pid.Voter{nissues}, space)
    for i in 1:nagents
    vi = pid.Voter(i, nissues, Tuple(1.,))
    pid.abm.add_agent_pos!(vi, model)
        end
    model
end

@testset "Am I generating the candidates correctly?" begin
  @test let m1 = simple_agents_model(10); ncandidates = 4
      pid.set_candidates!(ncandidates,m1)
      (pid.abm.allagents(m1) |>
          collect |>
          m-> sum(map(x->x.amIaCandidate, m))) == ncandidates
  end
end



function agents_model(nagents, nissues)
    space = abm.ContinuousSpace((1.,))
    properties = Dict(:current_incumbent => 0)
    model = abm.ABM(Voter{nissues}, space, properties = properties)
    for i in 1:nagents
    vi = Voter(i, nissues)
    abm.add_agent_pos!(vi, model)
    end
    return(model)
end


let m1 = agents_model(100); ncandidates = 4
      pid.set_candidates!(ncandidates,m1)
      (pid.abm.allagents(m1) |>
          collect |>
          m-> sum(map(x->x.amIaCandidate, m))) == ncandidates
  end


#= TODO: - Check if the id is actually a candidate, and Check if the candidate is actually the closest one
=#

# let model = agents_model(10)
#     testid = 4
#     pid.set_candidates!(2,model)
#     setfield!(model[testid], :pos, (0.5,))
#     for i in (pid.abm.allagents(model) |> collect) println(i) end
#     pid.get_closest_candidate(testid,model)
# end


# TODO: test if I am indeed choosing the most voted!

let model = agents_model(10)
    pid.set_candidates!(3,model)
    closest_candidates = Array{Int}(undef, pid.abm.nagents(model))
    for (fooindex,id) in enumerate(pid.abm.allids(model))
        closest_candidates[fooindex] = pid.get_closest_candidate(id,model)
    end
   pid.argmax(pid.proportionmap(closest_candidates)) == pid.getmostvoted(model)
end



# TODO: report there is something very wrong with the edistance proc in Agents.jl

#= for i in 1:100
let model = agents_model(1000)
    pid.set_candidates!(3,model)
    closest_candidates = Array{Int}(undef, pid.abm.nagents(model))
    somecandidate = pid.get_closest_candidate(1, model)
    model[somecandidate].pos = (0.5,)

    for (fooindex,id) in enumerate(pid.abm.allids(model))
        closest_candidates[fooindex] = pid.get_closest_candidate(id,model)
    end
   #println(model[somecandidate])
   #println(pid.proportionmap(closest_candidates))
    if pid.getmostvoted(model) != somecandidate

        println(model[pid.getmostvoted(model)], model[somecandidate])
        println(pid.proportionmap(closest_candidates))
        println(" ")
        #for i in collect(pid.abm.allagents(model)) println(i) end
    end
    #=So, this shows that it is not always true that the center will win
    #in this simulation because voters might be unevenly distributed!
    # TODO: Check if I am indeed picking the one that is closest to most agents
    #
=#
end
end
=#

let model = agents_model(10000)
    @time pid.set_candidates!(3,model)
    @time pid.getmostvoted(model)

    end


let model = agents_model(10000)
    @code_warntype pid.set_candidates!(3,model)
    @code_warntype pid.getmostvoted(model)

    end

# here is the last initilization step!
let model = agents_model(10000)
    pid.set_candidates!(3,model)
    new_incumbent = pid.getmostvoted(model)
    model.properties[:incumbent]= new_incumbent
    println(model.properties)
end

m1 = pid.initialize_model(1000, 3,2 )
@test m1.properties[:incumbent] == pid.getmostvoted(m1)
