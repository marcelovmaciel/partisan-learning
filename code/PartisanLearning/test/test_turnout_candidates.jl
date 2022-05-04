import Pkg

Pkg.activate("../PartisanLearning")

using Revise
import PartisanLearning as pl
import  Base.Filesystem as fl
using AlgebraOfGraphics
using CSV

# include("../test/visualize_model.jl")
# * Quick Vis 

m_params = pl.ModelParams(κ = 23, 
voter_pos_initializor = pl.hold(pl.sample_overlapping_2d_to_1d_hack,
pl.overlap_20_poss))

m = pl.initialize_model(m_params)

# ** Vis functions
fig,f = pl.single_interactive_vis(m); fig  



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


poodf = filter(:iter => x-> x==1.00, foodf)

filter(:step => x -> x==i, foodf)


myscatter(poodf)

myden() =   (data(foodf) * mapping(:f_ideal_point, :loyalty,  col = :sstep) * density(npoints = 50)) |> draw;

myscatter(df) =   (data(df) * mapping(:step, 
:loyalty, 
color = :f_ideal_point, 
size = :f_ideal_point => x-> x* 10) * (smooth() + visual(AlgebraOfGraphics.Scatter, colormap = :thermal))) |> draw;

mysurf() = data(poodf) *  mapping( :loyalty, :f_ideal_point, layout = :sstep) * AlgebraOfGraphics.density(npoints=50) *  visual(AlgebraOfGraphics.Surface) 

drawsurf() =  draw(mysurf(), axis=(type=AlgebraOfGraphics.Axis3, zticks=0:0.1:0.2, limits=(nothing, nothing, (0, 0.2))))

