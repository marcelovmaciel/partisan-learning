datapath = "../../../data"



function get_data_initial_dist(params, nsteps= 20)

    m = initialize_model(params)

    foursteps!(m)
    function loyalty(i) 
    loyalty(i.id,m)   # FIXME: I don't know why this is not working
    end 
    adata = [ f_ideal_point]

    mdata = [prop_pswitch,
             mean_loyalty,
             prop_crossvoting,
             get_2candidates_distance, 
             get_mean_contestant_eccentricity,
             get_incumbent_eccentricity]

    data_a,data_m= abm.run!(m,
                    abm.dummystep,
                    model_step!, nsteps;
                    adata,
                    mdata)
    return(data_a, data_m)
end


function collect_runs(params)
    holdermdf = DF.DataFrame([Int64[],
                          Float64[],
                          Int64[],
                          Float64[],
                          Float64[],
                          Float64[], 
                          Float64[],
                          Float64[]],
                         [:step,
                          :prop_pswitch,
                          :iter,
                          :mean_loyalty,
                          :prop_crossvoting,
                          :get_2candidates_distance, 
                          :get_mean_contestant_eccentricity,
                          :get_incumbent_eccentricity])
    holder_adf = DF.DataFrame(
        [Int64[], Int64[], Int64[], Float64[]],
        [:step, :id, :iter,  :f_ideal_point])
    
    @showprogress 1 "Running "  for i in 1:100
        adf,mdf = get_data_initial_dist(params)
        
        mdf[!,:iter] = [i for _ in 1:DF.nrow(mdf)]

        holdermdf = vcat(holdermdf,
                        mdf)
        adf[!,:iter] = [i for _ in 1:DF.nrow(adf)]
        
        println(names(holder_adf))
        println(names(adf))
        holder_adf = vcat(holder_adf,
                        adf)

    end

    return(holderadf, holdermdf)
end


function collect_per_overlap(whichdist, sampler)

ModelParams(voter_pos_initializor = hold(sampler, whichdist)) |>
        collect_runs
end
















# * FIXME: all code below is rotten
# only use it for inspiration!
function saltellidict(varnames::Vector{String}, bounds)
        Dict("names" => varnames,
             "bounds" => PythonCall.PyList(bounds),
             "num_vars" => length(varnames))
    end


function boundsdict_toparamsdf(varnames, bounds;samplesize= 2^8)
    saltelli = pyimport("SALib.sample.saltelli")

    boundsdict = saltellidict(varnames, bounds)

    problem_array = pyconvert(Array,
                              saltelli.sample(boundsdict, samplesize,
                                              calc_second_order = true ))
    foodf = DF.DataFrame()

    for (index,value) in enumerate(boundsdict["names"])
        setproperty!(foodf,Symbol(value), problem_array[:, index])
    end
    return(foodf)
end

function run_analysis_onRow(sim_cons,
                            row,
                            rowIdx,
                            whichIter,
                            threadId)
    m = initialize_model(sim_cons[:nagents], sim_cons[:nissues], row.ncandidates,
                             κ=  row.κ, δ =  row.δ)

    adata = [(a->(HaveIVotedAgainstMyParty(a,m)), x-> count(x)/m.properties[:nagents]),
             (a->(get_distance_IvsPartyCandidate(a,m)), d -> get_representativeness(d,m))]

    mdata = [normalized_ENP,
           x->x.properties[:incumbent_streak_counter].longest_streak[:streak_value],
             x-> x.properties[:party_switches][end]/x.properties[:nagents],
             get_incumbent_eccentricity]

    ad,md = abm.run!(m, abm.dummystep, model_step!, sim_cons[:niterations] ;
                       adata= adata,
                       mdata=mdata)#, when = collect_steps, when_model= collect_steps)

    alabels = ["step","Prop¬PartyId", "Representativeness"]
    mlabels = ["step", "NENP", "LongestIStreak", "PartySwitches", "Eccentricity"]
    DF.rename!(ad, alabels)
    DF.rename!(md,mlabels)
    df = hcat(ad, md, makeunique=true)
    DF.select!(df, DF.Not(:step_1))

    CSV.write(joinpath(datapath,
                       "row$(rowIdx)_iter$(whichIter)_thread$threadId.csv"),
              df)
    m = nothing
    df = nothing
end


function run_analysis(sim_cons, dm, threadId)

    ProgressMeter.@showprogress 5 "Simulating..."  for iter in 1:10
         for (rowidx,
                                             rowval) in enumerate(eachrow(dm))
            run_analysis_onRow(sim_cons,
                               rowval,
                               rowidx,
                               iter,
                               threadId)
        end
    end
end

get_mean_col_val(df,col,niterations)= StatsBase.mean(df[end-(niterations-1):end,col])

function get_system_measures(df)
    vars =  names(DF.select(df, DF.Not([:step])))

    Dict(Pair(var, get_mean_col_val(df, var, 20)) for var in vars)

end

readdf(dfname) = CSV.read(joinpath(datapath,dfname),
                          DF.DataFrame)

function system_measure_AtRepetions(whichParametization)
    data = readdir(datapath)

    map(get_system_measures∘readdf,
        filter(x-> ("row$(whichParametization)" == split(x, "_")[1]),
               data))
end

function get_ParametizationMeasuresMeans(whichparametrization)
    repetitionsvalues= DictionariesToDataFrame(system_measure_AtRepetions(whichparametrization))
    return(DF.mapcols(StatsBase.mean, repetitionsvalues))
end
