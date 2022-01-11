using InteractiveDynamics


function visualize_model(m)


  params = Dict(:κ => 0.0:1.:7, :switch => [:random, :plurality, :runoff])

  agent_colors(a) = a.id == m.properties[:incumbent]  ? :yellow : (a.amIaCandidate  ?  "#bf2642"  : "#2b2b33")
  agent_size(a) = a.id == m.properties[:incumbent]  ? 20 : (a.amIaCandidate ? 15 : 5)


  #higher kappa means agents will tend to vote more for
  #their partyid. Conversely, lower kappa means people vote more
  # for candidates closer to then rather than closer to their party
  # higher δ means parties sample further from their location.
  # both depended upon the underlying boundaries! think about that !!!

  adata = [#(a->(pla.HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
           (a->(pla.get_distance_IvsPartyCandidate(a,m)), d -> pla.get_representativeness(d,m))]

  mdata = [pla.normalized_ENP,
           x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
           x-> x.properties[:party_switches][end]/x.properties[:nagents],
           pla.get_incumbent_eccentricity]

  alabels = ["Rep"]
  mlabels = ["NENP", "IStreaks", "PSwitches", "ecc"]

  fig,adf,mdf = abm_data_exploration(m,
                                     pla.abm.dummystep,
                                     pla.model_step!,
                                     params;
                                     adata, mdata,
                                     alabels,mlabels,
                                     ac = agent_colors,
                                     as = agent_size, spu = 1 )


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
