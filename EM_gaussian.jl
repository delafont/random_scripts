
using MixtureModels     # fit_fmm
using Distributions     # Normal
using Gadfly
#using Winston           # (graphics)

obs = vcat( rand(Normal(2,1),1000), rand(Normal(6,1.5),400) )
r = fit_fmm(Normal{}, obs, 2, fmm_em(;display=:iter))

# r fields: (:mixture,:Q,:L,:niters,:converged,:objective)
r.mixture             # Mixture object, fields: (:K,:components,:pi)
r.mixture.components  # array of Normal objects

x = linspace(minimum(obs),maximum(obs),1000)


# Winston graphics
p = FramedPlot()
add(p, Points(obs,zeros(length(obs)), "type", "dot", "color", "blue"))
for comp in r.mixture.components
    y = pdf(comp,x)
    add(p, Points(x,y, "type", "dot", "color", "blue"))
end
Winston.display(p)


# Gadfly graphics
pl = plot(x=obs, Geom.point)
draw(PNG("myplot.svg"), pl)
