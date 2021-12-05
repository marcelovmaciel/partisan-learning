import Pkg
Pkg.activate("../")

using  Revise
using RParamMwe
using Distances

r = 1
for i in 1:20
    m = model_initialize(1000,2, 100, r )
    maximum([Distances.euclidean(m[cid].pos,m[nid].pos)
             for (cid,nid) in m.properties[:pointsdict]]) |> println
end



r = 2
for i in 1:20
    m = model_initialize(1000,2, 100, r )
    maximum([Distances.euclidean(m[cid].pos,m[nid].pos)
             for (cid,nid) in m.properties[:pointsdict]]) |> println
end


m.properties
