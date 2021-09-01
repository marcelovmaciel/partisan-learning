using Test
import PartisanLearning as pl
const is = pl.IssueSalience
const pid = pl.PartyId


 v1 = let amIaCandidate = false
    nissues = 1
    #pos = (rand(pid.distri.Uniform(0,1)),)
    pos = Tuple(1.)
    myPartyId = pos
    return(pid.Voter{nissues}(1,pos,amIaCandidate, myPartyId))
end
@test v1 == pid.Voter(1,1,Tuple(1.))

v2 = let space = pid.abm.ContinuousSpace((2.,))
    nissues = 1
    model = pid.abm.ABM(pid.Voter{nissues}, space)
    #pid.abm.add_agent!(Tuple(1.), model,false, (1.,))
    pid.abm.add_agent_pos!(pid.Voter(1, nissues, Tuple(1.,)), model)
    pid.abm.allagents(model) |> first
end

v3 = pid.Voter(1,1,Tuple(1.,))

function getstructvalues(S)
    map(s -> getfield(S,s), fieldnames(typeof(S)))
end

getstructvalues(v2) == getstructvalues(v3)

v4 = let space = pid.abm.ContinuousSpace((2.,))
    nissues = 1
    model = pid.abm.ABM(pid.Voter{nissues}, space)
    pid.abm.add_agent!(Tuple(1.), model,false, (1.,))
    pid.abm.allagents(model) |> first
end

getstructvalues(v4) == getstructvalues(v2)

#= So, I do not generally want to use the add_agent! procedure, because i tend
to control my agents' constructors =#





v2 == v3
