module RParamMwe

using Agents, Distances, Distributions

mutable struct  DummyAgent{n} <: AbstractAgent
    id::Int
    pos::NTuple{n,Float64}
end

function sample_Centerpoints_to_nearbydummies(npoints, model)
    ids = collect(allids(model))

    centerids= sample(ids,npoints,replace = false)

    centerpoints_nearbydummies = Dict(Pair(k,  -1) for k in centerids)

    return(centerpoints_nearbydummies)

end

function set_nearbies!(model)
    pointsdict = model.properties[:pointsdict]
    r = model.properties[:r]
    get_nearbypointid(pointid)= sample(collect(
        nearby_ids(model[pointid],
                       model, r,
                       exact = true)))

    centerpoint_nearbypoint_pairs = Dict(Pair(pid,
                                              get_nearbypointid(pid))
                                         for pid in keys(pointsdict))

    for (cid, nid) in centerpoint_nearbypoint_pairs
        pointsdict[cid] = nid
    end
end


function model_initialize(nagents,n, ncs,  r = 1)
    space = ContinuousSpace(ntuple(x -> float(10),n), periodic=false)

    model = ABM(DummyAgent{n}, space, properties = Dict(:n=> n,
                                                        :r => r,
                                                        :pointsdict => Dict{}()))

        for i in 1:nagents
        pos = Tuple(rand(Uniform(0,10), n))
                add_agent_pos!(DummyAgent{n}(i, pos), model)
                end


    model.properties[:pointsdict] = sample_Centerpoints_to_nearbydummies(ncs, model)
    set_nearbies!(model)
    return model
end

export DummyAgent, model_initialize, set_nearbies!
end # module
