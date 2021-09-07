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


#
# BUG: report there is something very wrong with the edistance proc in Agents.jl

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


@test let m1 = pid.initialize_model(1000, 3,2)
    m1.properties[:incumbent] == pid.getmostvoted(m1)
end


# TODO: ask on zulip how to check that without allocating
# @time reduce((x,y)-> x + y.amIaCandidate, pid.abm.allagents(m1), init = 0)

@testset "Check if reset_candidates keeps the incumbent as candidate" begin
    ncandidates = 3
    nissues = 1
    @test let m = pid.initialize_model(100,
                                       nissues,ncandidates)
        pid.reset_candidates!(m)
        reduce((x,y)-> x + y.amIaCandidate,
               pid.abm.allagents(m), init = 0) == 1
    end
    @test let m = pid.initialize_model(100,
                                       nissues,ncandidates)
        pid.reset_candidates!(m)
        pid.set_candidates!(ncandidates-1, m)
        reduce((x,y)-> x + y.amIaCandidate,
               pid.abm.allagents(m), init = 0) == 3
    end
end

ncandidates = 3
nissues = 1
m = pid.initialize_model(100,nissues, ncandidates)
pid.candidates_iteration_setup!(m)
# TODO check if candidates_iteration_setup! keeps the incumbent and add ncandidates - 1 other candidates



@testset "test if the candidate closest to agents differs from candidate closest to their partyid" begin
    @test let ncandidates = 10
        nissues = 1
        m = pid.initialize_model(1000,nissues, ncandidates)
        pid.candidates_iteration_setup!(m)
        any([begin
                 closest_to_me = pid.get_closest_candidate(agentid,m)
                 closest_to_myPartyId = pid.get_closest_candidate(agentid,
                                                                  m,
                                                                  :myPartyId)
            closest_to_me != closest_to_myPartyId
             end
             for agentid in pid.abm.allids(m)])
   end
end


@testset "test if the candidate an agent will vote for differs from candidate closest to their partyid" begin
    function  get_whoAgentVotesfor(agentid, model)
        κ = model.properties[:params].κ

        closest_to_me = get_closest_candidate(agentid,model)
        closest_to_myPartyId = get_closest_candidate(agentid, model, :myPartyId)
        whoillvotefor = closest_to_myPartyId
        two_candidates_distance = dist.euclidean(model[closest_to_me].pos,
                                                 model[closest_to_myPartyId].pos)
        if two_candidates_distance > κ
            whoillvotefor = closest_to_me
        end
        whoillvotefor != closest_to_myPartyId
    end
    @test let ncandidates = 10
        nissues = 1
        m = pid.initialize_model(1000,nissues, ncandidates, 0.1)
        pid.candidates_iteration_setup!(m)
        any([get_whoAgentVotesfor(agentid,m)
             for agentid in pid.abm.allids(m)])
    end
        @test let ncandidates = 2
            nissues = 1
            κ = 0.0001
            m = pid.initialize_model(1000,nissues, ncandidates, κ)
            pid.candidates_iteration_setup!(m)
            any([get_whoAgentVotesfor(agentid,m)
                 for agentid in pid.abm.allids(m)])
    end
end


# TODO: Write that as code!: IT MIGHT HAPPEN that people are very far from their partyid candidate!.

let ncandidates = 3
            nissues = 1
            κ = 0.5
            m = pid.initialize_model(1000,nissues, ncandidates, κ)
            pid.candidates_iteration_setup!(m)
            any([pid.get_whoAgentVotesfor(agentid,m)
                 for agentid in pid.abm.allids(m)])
end

# this shows that this procedure might have a huge impact !
let ncandidates = 3
    nissues = 1
    κ = 0.1
    m = pid.initialize_model(1000,nissues, ncandidates, κ)
    pid.candidates_iteration_setup!(m)
    println(m[1].myPartyId)
    pid.update_partyid!(1,m)
    println(m[1].myPartyId)

end


@test let ncandidates = 3
    nissues = 3
    κ = 0.1
    m = pid.initialize_model(1000,nissues, ncandidates, κ)
    pid.candidates_iteration_setup!(m)
     for i in pid.abm.allids(m)
        pid.update_partyid!(i,m)
     end

    agentid = 1
    ρ = 0.1
    length(pid.get_median_neighborhoodPid(agentid,m, ρ)) == nissues
end



# FIXME: Investigate the performance of model_step!
# let ncandidates = 3
#     nissues = 3
#     κ = 0.1
#     m = pid.initialize_model(1000,nissues, ncandidates, κ)
#     @time pid.model_step!(m)

#      # for i in pid.abm.allids(m)
#      #    pid.update_partyid!2(i,m)
#      # endd
# end
