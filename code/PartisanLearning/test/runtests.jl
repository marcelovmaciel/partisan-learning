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


agents_model(nagents) = let space = pid.abm.ContinuousSpace((2.,))
    nissues = 1
    model = pid.abm.ABM(pid.Voter{nissues}, space)
    for i in 1:nagents
    vi = pid.Voter(i, nissues)
    pid.abm.add_agent_pos!(vi, model)
        end
    model
end


@testset "Am I getting the closest candidate as designed? Or is there a bug" begin
    @test let model = agents_model(10)
        testid = 4
        pid.set_candidates!(2,model)
        setfield!(model[testid], :pos, (0.5,))
        for i in (pid.abm.allagents(model) |> collect) println(i) end
        dummypos = 10.; dummyid = -2
        for i in
            pid.abm.allids(model)     #pid.abm.nearby_ids(model[testid].pos, model)
            if !model[i].amIaCandidate
                continue
            else
                if pid.abm.edistance(testid, i, model) < dummypos
                    dummypos = pid.abm.edistance(testid, i, model)
                    dummyid = i
                end
            end
        end
        dummyid == pid.get_closest_candidate(testid, model)
    end
end
# let model = agents_model(10)
#     testid = 4
#     pid.set_candidates!(2,model)
#     setfield!(model[testid], :pos, (0.5,))
#     for i in (pid.abm.allagents(model) |> collect) println(i) end
#     pid.get_closest_candidate(testid,model)
# end



#= TODO:
- Check if the id is actually a candidate
- Check if the candidate is actually the closest one
=#
