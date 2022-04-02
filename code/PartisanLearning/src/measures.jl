
# ** Data Collection
#=

FIXME: Test what happens with the  proportion of voters who voted against PartyId candidate

Maybe also some measures of the distribution? Who knows....

• for how long candidate remains winning.
• The proportion of voters voted for someone who was not their party
candidate.
• What is the distance between the candidate the agent voted for and
the candidate of their party.
=#


#adata = [(a->(HaveIVotedAgainstMyParty(a,m)), +)]

function get_distance_MyCandidatevsPartyCandidate(agentid, model)


    closest_to_me = get_closest_candidate(agentid,model)[1]
    mypartycandidate = model.properties[:parties_candidateid_ppos_δ][model[agentid].myPartyId][:partycandidate]
    two_candidates_distance = dist.euclidean(model[closest_to_me].pos,
                                             model[mypartycandidate].pos)
    return(two_candidates_distance)

end

function get_distance_IvsPartyCandidate(agentid,model)
    mypartycandidate = model.properties[:parties_candidateid_ppos_δ][model[agentid].myPartyId][:partycandidate]
    return(dist.euclidean(model[agentid].pos,
                          model[mypartycandidate].pos))
end


get_distance_IvsPartyCandidate(agent::Voter,model) = get_distance_IvsPartyCandidate(agent.id,model)
get_distance_MyCandidatevsPartyCandidate(agent::Voter, model) = get_distance_MyCandidatevsPartyCandidate(agent.id, model)


function get_representativeness(m::abm.ABM)
    -sum([get_distance_IvsPartyCandidate(i,m) for i in abm.allids(m)])/m.properties[:nagents]
end

function get_representativeness(v,m)
    -sum(v)/m.properties[:nagents]
end

function get_partyshare(m)
    proportionmap([m.properties[:voterBallotTracker][agentid][end]
                   for agentid in abm.allids(m)])
end

function get_ENP(m)
    shares = get_partyshare(m)
    return(1/sum([i^2 for i in collect(values(shares))]))
end


function normalized_ENP(m)
    get_ENP(m)/m.properties[:ncandidates]
end


# x->x[x.properties[:incumbent_party]].myPartyId, # get incumbent_party

function get_incumbent_eccentricity(m)
      dist.euclidean(m[m.properties[:incumbent_party]].pos,
                   m.properties[:median_pos])
end


function get_mean_contestant_eccentricity(m)
    contestants = filter(pid-> pid != m.properties[:incumbent_party],
                         m.properties[:parties_ids])
    (contestants .|>
      (contestant -> dist.euclidean(m[contestant].pos,
                                    m.properties[:median_pos])) |>
                                        mean)
end



function get_party_supporters_mean(pid, m, whichdim = 1 )
    party_supporters = get_parties_supporters(m)[pid]
    mean(vcat(map(x->collect(m[x].pos[whichdim]), party_supporters)...))
end




function get_2candidates_distance(m)
    parties = m.properties[:parties_ids]
    candidate1_pos = m[m.properties[:parties_candidateid_ppos_δ][parties[1]][:partycandidate]].pos
    candidate2_pos = m[m.properties[:parties_candidateid_ppos_δ][parties[2]][:partycandidate]].pos
    dist.euclidean(candidate1_pos,
                   candidate2_pos)
end


function mean_loyalty(m)
[get_keep_party_id_prob(i,m) for i in abm.allids(m)] |>  mean
        end

loyalty(i) = get_keep_party_id_prob(i.id,m)
f_ideal_point(i) = i.pos[1]

function prop_pswitch(m)
    m.properties[:party_switches][end]/m.properties[:nagents]
end

function prop_crossvoting(m)
    m.properties[:cross_voting][end]/m.properties[:nagents]
end
