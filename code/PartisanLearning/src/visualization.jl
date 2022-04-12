# TODO: one function should dispatch interactive exploration
# Another function should dispatch sweeping exploration

# * One type of model interactive vis

agent_size(a,m) = (a.id in m.properties[:parties_ids] ?  37 : (a.amIaCandidate ? 35 : 5))
agent_marker(a,m) = if a.id in m.properties[:parties_ids] '♠' else '∘' end

function agent_colors_2parties(a,m)

    p1 = m.properties[:parties_ids][1]
    p2 = m.properties[:parties_ids][2]

    if a.id == m.properties[:incumbent_party]
        :yellow
    elseif a.myPartyId == p1
        :red
    elseif a.myPartyId == p2
        :blue
    end
end

# params = Dict(:κ => 0.0:1:100.0,
    #           :switch => [:random, :plurality, :runoff],
    #           :δ=> 0.0:5:70,
    #           :ω => 0.0:0.1:1.0,
    #           :kappa_switch => [:off, :on]) FIXME: this breaks the code

function single_interactive_vis(m)
    foursteps!(m)
    az(a) = agent_size(a,m)
    am(a) = agent_marker(a,m)
    ac(a) = agent_colors_2parties(a,m)
    println("got here")
    adata = [(a->(get_distance_IvsPartyCandidate(a,m)), d -> get_representativeness(d,m))]

    mdata =  [get_2candidates_distance, prop_pswitch, mean_loyalty]

    alabels = ["rep"]
    mlabels = ["dist c1c2", "prop_switch","mean_loyalty"]
    fig,pl2 = abmexploration(m;
                             agent_step! = abm.dummystep,
                             model_step! = model_step!,
                            # adata, # BUG: this breaks the code !
                              mdata,
                             # alabels,
                              mlabels,
                             ac = ac,
                             as = az,
                             am = am,
                             spu = 1:100)

    return(fig, pl2)
end