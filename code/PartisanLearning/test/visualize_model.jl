using InteractiveDynamics


function visualize_model(m)
    # FIXME: getting a problem now that initial condition is only initial lol
    pla.model_step!(m)

    params = Dict(:κ => 0.0:10:70, :switch => [:random, :plurality, :runoff],
                  :δ => 0.0:10:70,
                  :ω => 0.0:0.1:1.0,
                  :kappa_switch => [:on, :off])

    agent_colors(a) = a.id == m.properties[:incumbent_party]  ? :yellow : (a.amIaCandidate  ?  "#bf2642"  : "#2b2b33")
    agent_size(a) = a.id == m.properties[:incumbent_party]  ? 25 : (a.id in m.properties[:parties_ids] ? 20 : (a.amIaCandidate ? 15 : 5))
    agent_marker(a) = if a.id in m.properties[:parties_ids] '♠' else '∘' end

  #higher kappa means agents will tend to vote more for
  #their partyid. Conversely, lower kappa means people vote more
  # for candidates closer to then rather than closer to their party
  # higher δ means parties sample further from their location.
  # both depended upon the underlying boundaries! think about that !!!

  adata = [#(a->(pla.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
           (a->(pla.get_distance_IvsPartyCandidate(a,m)), d -> pla.get_representativeness(d,m))]

  mdata = [
           x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
           x-> x.properties[:party_switches][end]/x.properties[:nagents],
           x-> x.properties[:incumbent_streak_counter].has_switchedlist[end],
           pla.get_incumbent_eccentricity,
           pla.get_mean_contestant_eccentricity]
#pla.normalized_ENP,
  alabels = ["Rep"]
  mlabels = [ "IStreaks", "PSwitches","iSwitch", "ie", "ce"]
#"NENP",
  fig,adf,mdf = abm_data_exploration(m,
                                     pla.abm.dummystep,
                                     pla.model_step!,
                                     params;
                                     adata, mdata,
                                     alabels,mlabels,
                                     ac = agent_colors,
                                     as = agent_size,
                                     am = agent_marker,  spu = 1 )


  partyproportions = pla.proportionmap([v for (k,v) in m.properties[:voters_partyids]])
  partyids = collect(keys(partyproportions))
  foo = Observable([partyproportions[p] for p in partyids])


  function newstep(m, foo = foo)
      pla.model_step!(m)
      partyproportions = pla.proportionmap([v for (k,v) in m.properties[:voters_partyids]])
      foo[] = [partyproportions[p] for p in partyids]
  end

  fig,adf,mdf = abm_data_exploration(m,
                                     pla.abm.dummystep,
                                     newstep,
                                     params;
                                     adata, mdata,
                                     alabels,mlabels,
                                     ac = agent_colors,
                                     as = agent_size, spu = 1 )

  barplot(fig[1,3], 1:length(partyids), foo, title = "Proportion Between Parties" )


end


haveiswitched(agent::pla.Voter, m) = haveiswitched(agent.id,m)

function visualize_simpler(m)
    # FIXME: getting a problem now that initial condition is only initial lol
    pla.model_step!(m)

    params = Dict(:κ => 0.0:5:100.0, :switch => [:random, :plurality, :runoff],
                  :δ => 0.0:10:70,
                  :ω => 0.0:0.1:1.0,
                  :kappa_switch => [:off, :on])

    function agent_colors(a)

        if a.id == m.properties[:incumbent_party]
            :yellow
        elseif a.myPartyId == m.properties[:parties_ids][1]
            :red
        else :blue
        end

        end
    agent_size(a) = a.id == m.properties[:incumbent_party]  ? 25 : (a.id in m.properties[:parties_ids] ? 20 : (a.amIaCandidate ? 15 : 5))
    agent_marker(a) = if a.id in m.properties[:parties_ids] '♠' else '∘' end

  #higher kappa means agents will tend to vote more for
  #their partyid. Conversely, lower kappa means people vote more
  # for candidates closer to then rather than closer to their party
  # higher δ means parties sample further from their location.
  # both depended upon the underlying boundaries! think about that !!!

  # adata = [#(a->(pla.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
  #          (a->(pla.get_distance_IvsPartyCandidate(a,m)), d -> pla.get_representativeness(d,m))]


    adata = [(i->pla.get_keep_party_id_prob(i.id,m), pla.mean)]
    mdata = [ x-> x.properties[:party_switches][end]/x.properties[:nagents]]
            #pla.get_incumbent_eccentricity,
            #pla.get_mean_contestant_eccentricity]

  #          x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
  #
  #          x-> x.properties[:incumbent_streak_counter].has_switchedlist[end],
  #          pla.get_incumbent_eccentricity,
  #          pla.get_mean_contestant_eccentricity]
#pla.normalized_ENP,
    alabels = ["keep party prob"]
    mlabels = ["PSwitches"]
  # mlabels = [ "IStreaks", "PSwitches","iSwitch", "ie", "ce"]
#"NENP",
  fig,adf,mdf = abm_data_exploration(m,
                                     pla.abm.dummystep,
                                     pla.model_step!,
                                     params;
                                     adata,
                                     mdata,
                                     alabels,
                                     mlabels,
                                     ac = agent_colors,
                                     as = agent_size,
                                     am = agent_marker,  spu = 1 )


  #keep_pid_probs = [pla.get_keep_party_id_prob(i,m) for i in  allids(m)]




  # foo = Observable(keep_pid_probs)
  # n = Observable([0.])



  # function newstep(m, foo = foo)
  #     pla.model_step!(m)
  #     keep_pid_probs = [pla.get_keep_party_id_prob(i,m) for i in  allids(m)]
  #     foo[] = keep_pid_probs

  # end

  # fig,adf,mdf = abm_data_exploration(m,
  #                                    pla.abm.dummystep,
  #                                    newstep,
  #                                    params;
  #                                       ac = agent_colors,
  #                                    as = agent_size, spu = 1 )
  #   hist(fig[1,3], foo)



end






function visualize_simpler_pluscustom(m)
    # FIXME: getting a problem now that initial condition is only initial lol
    pla.model_step!(m)

    params = Dict(:κ => 0.0:5:100.0, :switch => [:random, :plurality, :runoff],
                  :δ => 0.0:10:70,
                  :ω => 0.0:0.1:1.0,
                  :kappa_switch => [:off, :on])

    function agent_colors(a)

        if a.id == m.properties[:incumbent_party]
            :yellow
        elseif a.myPartyId == m.properties[:parties_ids][1]
            :red
        else :blue
        end

        end
    agent_size(a) = a.id == m.properties[:incumbent_party]  ? 25 : (a.id in m.properties[:parties_ids] ? 20 : (a.amIaCandidate ? 15 : 5))
    agent_marker(a) = if a.id in m.properties[:parties_ids] '♠' else '∘' end

  #higher kappa means agents will tend to vote more for
  #their partyid. Conversely, lower kappa means people vote more
  # for candidates closer to then rather than closer to their party
  # higher δ means parties sample further from their location.
  # both depended upon the underlying boundaries! think about that !!!

  # adata = [#(a->(pla.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
  #          (a->(pla.get_distance_IvsPartyCandidate(a,m)), d -> pla.get_representativeness(d,m))]


    adata = [(i->pla.get_keep_party_id_prob(i.id,m), pla.mean)]
    mdata = [ x-> x.properties[:party_switches][end]/x.properties[:nagents]]

  #          x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
  #
  #          x-> x.properties[:incumbent_streak_counter].has_switchedlist[end],
  #          pla.get_incumbent_eccentricity,
  #          pla.get_mean_contestant_eccentricity]
#pla.normalized_ENP,
    alabels = ["keep party prob"]
    mlabels = ["PSwitches"]
  # mlabels = [ "IStreaks", "PSwitches","iSwitch", "ie", "ce"]
#"NENP",
  fig,adf,mdf = abm_data_exploration(m,
                                     pla.abm.dummystep,
                                     pla.model_step!,
                                     params;
                                     # adata,
                                     # mdata,
                                     # alabels,
                                     # mlabels,
                                     ac = agent_colors,
                                     as = agent_size,
                                     am = agent_marker,  spu = 1 )



    keep_pid_probs = [Observable([convert(Float32,pla.get_keep_party_id_prob(i,m))])
                      for i in  allids(m)]
    N = Observable([Float32(0.)])


  function newstep(m, keep_pid_probs = keep_pid_probs, N = N)
      pla.model_step!(m)
      probs = [convert(Float32,pla.get_keep_party_id_prob(i,m)) for i in  allids(m)]
      for (i,v) in enumerate(probs)
     keep_pid_probs[i][] = [v]
      end
      N[] = push!(N.val,N.val[end]+ 1.)
  end

  fig,adf,mdf = abm_data_exploration(m,
                                     pla.abm.dummystep,
                                     newstep,
                                     params;
                                        ac = agent_colors,
                                     as = agent_size, spu = 1 )


    scatter(fig[1,3], N, keep_pid_probs[1])
    lines!(fig[1,3],N, keep_pid_probs[1])
    for i in keep_pid_probs[2:end]
        scatter!(fig[1,3], N, i)
        lines!(fig[1,3], N, i)
    end
    autolimits!(fig[1,3])
return(keep_pid_probs, N)

end





function visualize_m(m)
    pla.model_step!(m)
    fig,adf,mdf = abm_data_exploration(m,
                                       pla.abm.dummystep,
                                       pla.model_step!)
    var_wanna_observe = [Observable([pla.get_keep_party_id_prob(i,m)])
                                    for i in  allids(m)]
    nsteps = Observable([0.])
    function newstep(m, var_wanna_observe = var_wanna_observe, nsteps = nsteps)
        pla.model_step!(m)
        var_values = [pla.get_keep_party_id_prob(i,m) for i in  allids(m)]
        for (i,v) in enumerate(var_values)
            push!(var_wanna_observe[i].val, v)
        end
        nsteps[] = push!(nsteps.val,nsteps.val[end]+ 1.)
    end
    fig,adf,mdf = abm_data_exploration(m, dummystep, newstep)
    ax = Axis(fig[1,2])
    for i in var_wanna_observe[1:end]
        lines!(ax, nsteps, i, linewidth = 0.5)
        autolimits!(ax)
    end
    autolimits!(ax)
end
