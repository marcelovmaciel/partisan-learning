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


# * One type of Model sweeping vis


function get_data_initial_dist(params, nsteps= 20)

    m = initialize_model(params)

    foursteps!(m)


    adata = [f_ideal_point]# ,
             # loyalty] FIXME: this breaks the code!

    mdata = [prop_pswitch,
             mean_loyalty,
             prop_crossvoting,
             get_2candidates_distance]

    _,data_m= run!(m,
                    abm.dummystep,
                    model_step!, nsteps;
                    adata,
                    mdata)
    return(data_m)
end



function collect_runs(params)

    holderdf = DataFrame([Int64[],
                          Float64[],
                          Int64[],
                          Float64[],
                          Float64[],
                          Float64[]],
                         [:step,
                          :prop_pswitch,
                          :iter,
                          :mean_loyalty,
                          :prop_crossvoting,
                          :get_2candidates_distance])

@showprogress 1 "Running "  for i in 1:100
        df = get_data_initial_dist(params)
        df[!,:iter] = [i for _ in 1:21]

        holderdf = vcat(holderdf,
                        df)
    end

    return(holderdf)
end


function collect_per_overlap(whichdist)
    (whichdist |>
        overlap_initializor |>
        x-> ModelParams(voter_pos_initializor = x) |>
        collect_runs)
end
