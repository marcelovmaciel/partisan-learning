import Pkg

Pkg.activate("../")

using Revise
import PartisanLearning as pl

# include("../test/visualize_model.jl")

m_params = pl.ModelParams()

m = pl.initialize_model(m_params)

# * Vis functions
fig,f = single_interactive_vis(m)


overlap20df = collect_per_overlap(pl.overlap_20_poss)

overlap50df = collect_per_overlap(pl.overlap_50_poss)

overlap80df = collect_per_overlap(pl.overlap_80_poss)

function concat_overlap_data(f,s,t)

    for (i,j) in zip([f,s,t], [:o20, :o50, :o80])
        println(j)
        @show first(i)

        df_size =  (i |> size)[1]
        i[!,:type] =[j for _ in 1:df_size]
    end
    gdf = vcat(f,s,t)
    return(gdf)
end


gdf =  concat_overlap_data(overlap20df,
                           overlap50df,
                           overlap80df)

plt = data(gdf) * mapping(:step,
                          :get_2candidates_distance,
                          color = :type => nonnumeric,
                          dodge = :type => nonnumeric) * visual(BoxPlot);

foo = draw(plt)

stringtosave = "./plots/overlapping_2cdist_kappa1.png"

AlgebraOfGraphics.save(stringtosave, foo, px_per_unit = 5)
