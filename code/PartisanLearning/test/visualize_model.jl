using InteractiveDynamics


 function visualize_model(m)
    # FIXME: getting a problem now that initial condition is only initial lol


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

  adata = [#(a->(pl.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
           (a->(pl.get_distance_IvsPartyCandidate(a,m)), d -> pl.get_representativeness(d,m))]

  mdata = [
           x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
           x-> x.properties[:party_switches][end]/x.properties[:nagents],
           x-> x.properties[:incumbent_streak_counter].has_switchedlist[end],
           pl.get_incumbent_eccentricity,
           pl.get_mean_contestant_eccentricity]
#pl.normalized_ENP,
  alabels = ["Rep"]
  mlabels = [ "IStreaks", "PSwitches","iSwitch", "ie", "ce"]
#"NENP",
  fig,adf,mdf = abm_data_exploration(m,
                                     pl.abm.dummystep,
                                     pl.model_step!,
                                     params;
                                     adata, mdata,
                                     alabels,mlabels,
                                     ac = agent_colors,
                                     as = agent_size,
                                     am = agent_marker,  spu = 1 )


  partyproportions = pl.proportionmap([v for (k,v) in m.properties[:voters_partyids]])
  partyids = collect(keys(partyproportions))
  foo = Observable([partyproportions[p] for p in partyids])


  function newstep(m, foo = foo)
      pl.model_step!(m)
      partyproportions = pl.proportionmap([v for (k,v) in m.properties[:voters_partyids]])
      foo[] = [partyproportions[p] for p in partyids]
  end

  fig,adf,mdf = abm_data_exploration(m,
                                     pl.abm.dummystep,
                                     newstep,
                                     params;
                                     adata, mdata,
                                     alabels,mlabels,
                                     ac = agent_colors,
                                     as = agent_size, spu = 1 )

  barplot(fig[1,3], 1:length(partyids), foo, title = "Proportion Between Parties" )


end


haveiswitched(agent::pl.Voter, m) = haveiswitched(agent.id,m)


function visualize_simpler(m)
    # FIXME: getting a problem now that initial condition is only initial lol
    pl.model_step!(m)
     pl.model_step!(m)
     pl.model_step!(m)
    pl.model_step!(m)


    params = Dict(:κ => 0.0:1:100.0, :switch => [:random, :plurality, :runoff],
                  :δ => 0.0:5:70,
                  :ω => 0.0:0.1:1.0,
                  :kappa_switch => [:off, :on])

    p1 = m.properties[:parties_ids][1]
    p2 = m.properties[:parties_ids][2]

    function agent_colors(a)

        if a.id == m.properties[:incumbent_party]
            :yellow
        elseif a.myPartyId == p1
            :red
        elseif  a.myPartyId == p2
            :blue
        end

        end
    agent_size(a) = a.id == m.properties[:incumbent_party]  ? 25 : (a.id in m.properties[:parties_ids] ? 20 : (a.amIaCandidate ? 15 : 5))
    agent_marker(a) = if a.id in m.properties[:parties_ids] '♠' else '∘' end

  #higher kappa means agents will tend to vote more for
  #their partyid. Conversely, lower kappa means people vote more
  # for candidates closer to then rather than closer to their party
  # higher δ means parties sample further from their location.
  # both depended upon the underlying boundaries! think about that !!!

   adata = [(a->(pl.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents])]
  #          (a->(pl.get_distance_IvsPartyCandidate(a,m)), d -> pl.get_representativeness(d,m))]


  #  adata = [(i->pl.get_keep_party_id_prob(i.id,m), pl.mean)]
    # mdata = [m-> m[m.properties[:parties_ids][1] ].pos[1],
      #       m -> m[m.properties[:parties_ids][2] ].pos[1],
       #      pl.get_2candidates_distance]
        #[ x-> x.properties[:party_switches][end]/x.properties[:nagents]]
            #pl.get_incumbent_eccentricity,
            #pl.get_mean_contestant_eccentricity]

  #          x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
  #
  #          x-> x.properties[:incumbent_streak_counter].has_switchedlist[end],
  #          pl.get_incumbent_eccentricity,
  #          pl.get_mean_contestant_eccentricity]
#pl.normalized_ENP,
    alabels = ["prop voting against party"]
    # mlabels = ["Party1Pos", "Party2Pos", "dist(c1,c2)"]
  # mlabels = [ "IStreaks", "PSwitches","iSwitch", "ie", "ce"]
#"NENP",
  fig,adf,mdf = abm_data_exploration(m,
                                     pl.abm.dummystep,
                                     pl.model_step!,
                                     params;
                                     adata,
        #                             mdata,
                                     alabels,
      #                               mlabels,
                                     ac = agent_colors,
                                     as = agent_size,
                                     am = agent_marker,  spu = 1:100 )


  #keep_pid_probs = [pl.get_keep_party_id_prob(i,m) for i in  allids(m)]




  # foo = Observable(keep_pid_probs)
  # n = Observable([0.])



  # function newstep(m, foo = foo)
  #     pl.model_step!(m)
  #     keep_pid_probs = [pl.get_keep_party_id_prob(i,m) for i in  allids(m)]
  #     foo[] = keep_pid_probs

  # end

  # fig,adf,mdf = abm_data_exploration(m,
  #                                    pl.abm.dummystep,
  #                                    newstep,
  #                                    params;
  #                                       ac = agent_colors,
  #                                    as = agent_size, spu = 1 )
  #   hist(fig[1,3], foo)



end


function visualize_simpler_1d(m)
    # FIXME: getting a problem now that initial condition is only initial lol
    pl.model_step!(m)

    params = Dict(:κ => 0.0:1:100.0, :switch => [:random, :plurality, :runoff],
                  :δ => 0.0:5.:70,
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

  # adata = [#(a->(pl.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
  #          (a->(pl.get_distance_IvsPartyCandidate(a,m)), d -> pl.get_representativeness(d,m))]


    adata = [(i->pl.get_keep_party_id_prob(i.id,m), pl.mean)]
    mdata = [ x-> x.properties[:party_switches][end]/x.properties[:nagents]]
            #pl.get_incumbent_eccentricity,
            #pl.get_mean_contestant_eccentricity]

  #          x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
  #
  #          x-> x.properties[:incumbent_streak_counter].has_switchedlist[end],
  #          pl.get_incumbent_eccentricity,
  #          pl.get_mean_contestant_eccentricity]
#pl.normalized_ENP,
    alabels = ["keep party prob"]
    mlabels = ["PSwitches"]
  # mlabels = [ "IStreaks", "PSwitches","iSwitch", "ie", "ce"]
#"NENP",



  fig,adf,mdf = abm_data_exploration(m,
                                     pl.abm.dummystep,
                                     pl.model_step!,
                                     params;
                                     adata,
                                     mdata,
                                     alabels,
                                     mlabels,
                                     ac = agent_colors,
                                     as = agent_size,
                                     am = agent_marker,  spu = 1:100 )


  #keep_pid_probs = [pl.get_keep_party_id_prob(i,m) for i in  allids(m)]



    ax1 = Axis(fig[1,2])
    hideydecorations!(ax1, ticks = false)
    hidespines!(ax1)
  # foo = Observable(keep_pid_probs)
  # n = Observable([0.])



  function newstep(m, foo = foo)
      pl.model_step!(m)
      keep_pid_probs = [pl.get_keep_party_id_prob(i,m) for i in  allids(m)]
      foo[] = keep_pid_probs

  end

  # fig,adf,mdf = abm_data_exploration(m,
  #                                    pl.abm.dummystep,
  #                                    newstep,
  #                                    params;
  #                                       ac = agent_colors,
  #                                    as = agent_size, spu = 1 )
  #   hist(fig[1,3], foo)



end




function visualize_noslider(m, spu = 1:100)
    # FIXME: getting a problem now that initial condition is only initial lol
    pl.model_step!(m)

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

  # adata = [#(a->(pl.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
  #          (a->(pl.get_distance_IvsPartyCandidate(a,m)), d -> pl.get_representativeness(d,m))]


    adata = [(i->pl.get_keep_party_id_prob(i.id,m), pl.mean)]
    mdata = [ x-> x.properties[:party_switches][end]/x.properties[:nagents]]
            #pl.get_incumbent_eccentricity,
            #pl.get_mean_contestant_eccentricity]

  #          x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
  #
  #          x-> x.properties[:incumbent_streak_counter].has_switchedlist[end],
  #          pl.get_incumbent_eccentricity,
  #          pl.get_mean_contestant_eccentricity]
#pl.normalized_ENP,
    alabels = ["keep party prob"]
    mlabels = ["PSwitches"]
  # mlabels = [ "IStreaks", "PSwitches","iSwitch", "ie", "ce"]
#"NENP",
  fig,adf,mdf = abm_data_exploration(m,
                                     pl.abm.dummystep,
                                     pl.model_step!;
                                     adata,
                                     mdata,
                                     alabels,
                                     mlabels,
                                     ac = agent_colors,
                                     as = agent_size,
                                     am = agent_marker,  spu = spu )





end






function visualize_simpler_pluscustom(m)
    # FIXME: getting a problem now that initial condition is only initial lol
    pl.model_step!(m)

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

  # adata = [#(a->(pl.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
  #          (a->(pl.get_distance_IvsPartyCandidate(a,m)), d -> pl.get_representativeness(d,m))]


    adata = [(i->pl.get_keep_party_id_prob(i.id,m), pl.mean)]
    mdata = [ x-> x.properties[:party_switches][end]/x.properties[:nagents]]

  #          x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
  #
  #          x-> x.properties[:incumbent_streak_counter].has_switchedlist[end],
  #          pl.get_incumbent_eccentricity,
  #          pl.get_mean_contestant_eccentricity]
#pl.normalized_ENP,
    alabels = ["keep party prob"]
    mlabels = ["PSwitches"]
  # mlabels = [ "IStreaks", "PSwitches","iSwitch", "ie", "ce"]
#"NENP",
  fig,adf,mdf = abm_data_exploration(m,
                                     pl.abm.dummystep,
                                     pl.model_step!,
                                     params;
                                     # adata,
                                     # mdata,
                                     # alabels,
                                     # mlabels,
                                     ac = agent_colors,
                                     as = agent_size,
                                     am = agent_marker,  spu = 1 )



    keep_pid_probs = [Observable([convert(Float32,pl.get_keep_party_id_prob(i,m))])
                      for i in  allids(m)]
    N = Observable([Float32(0.)])


  function newstep(m, keep_pid_probs = keep_pid_probs, N = N)
      pl.model_step!(m)
      probs = [convert(Float32,pl.get_keep_party_id_prob(i,m)) for i in  allids(m)]
      for (i,v) in enumerate(probs)
     keep_pid_probs[i][] = [v]
      end
      N[] = push!(N.val,N.val[end]+ 1.)
  end

  fig,adf,mdf = abm_data_exploration(m,
                                     pl.abm.dummystep,
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
    pl.model_step!(m)
    fig,adf,mdf = abm_data_exploration(m,
                                       pl.abm.dummystep,
                                       pl.model_step!)
    var_wanna_observe = [Observable([pl.get_keep_party_id_prob(i,m)])
                                    for i in  allids(m)]
    nsteps = Observable([0.])
    function newstep(m, var_wanna_observe = var_wanna_observe, nsteps = nsteps)
        pl.model_step!(m)
        var_values = [pl.get_keep_party_id_prob(i,m) for i in  allids(m)]
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
