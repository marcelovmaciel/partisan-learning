using  Revise
using RParamMwe
using Distances


for i in 1:20
    m = model_initialize(1000,5, 100, 1 )
    maximum([Distances.euclidean(m[cid].pos,m[nid].pos)
             for (cid,nid) in m.properties[:pointsdict]]) |> println
end
