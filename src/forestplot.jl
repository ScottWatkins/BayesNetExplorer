"""
    forestplot(n::Vector, m::Vector, e1::Vector, e2::Vector)

where n contains names, m contains point estimates, e1 (& e2) contain the errors (m +/- e) values, and xlabel is the x-axis label. If error vectors are the actual CI plot values, e1 should have the lower values and e2 the upper values and set CI=true. otherwise use e1.

    Options
        markercolor        series1
        markercolor2       series2
        markersize         series1
        markersize2        series2
        minorticks         5
        scale_x            set x-scale        [:identity|:log10]
        xtickvals          set x-axis ticks   []
        xlims              set x limits       (1,16)
        xlabel             x-axis label       ""
        legend             legend placement   :best
        label1             series 1 label     ""
        label2             series 2 label     "" 
        CI                 error values       false
"""
function forestplot(n::Vector, m::Vector, e1::Vector, e2::Vector; xlabel::String="Odds Ratio", CI::Bool=false, markercolor::Symbol=:midnightblue, markershape::Symbol=:circle, markersize::Int64=4, m2::Vector=[], m2e1::Vector=[], m2e2::Vector=[], markercolor2::Symbol=:firebrick, markershape2::Symbol=:circle, markersize2::Int64=7, scale_x=:identity, legend=:best, xtickvals::Array=[], xlims::Tuple{Number,Number}=(0,10), minorticks=5, label1::String="", label2::String="")

    length(n) == length(m) == length(e1) == length(e2) ? println("Input vector length equal...pass.") : error("Vector length error.")
    println("Plotting...")

    if CI == true
        le = abs.(m .- e1)
        ue = abs.(e2 .- m)
    else
        if length(e2) > 0
            if e1 != e2
                println("Found values for e2 but setting e2 == e1. Set CI=true to use e2 values.")
            end
        end
        le = e1
        ue = e1
    end
    
    xlabel= xlabel * "\n"
    
    ypad = 0.72
    err = Tuple.(zip(le,ue))

    if length(m2) > 0
        y = collect(1.0:1.0:length(m))
        y1 = y .- 0.2
        y2 = y .+ 0.2
    else
        y = collect(1:1:length(m))
        y1 = y
    end

    plt = scatter(m, y1, xaxis=scale_x, xerr=err,
            #markershape=markershape, #can't use markershape due to bug.
            markersize = markersize,
            markeralpha = 0.99,
            markercolor = :black,
            markerstrokewidth = 1,
            markerstrokealpha = 0.99,
            markerstrokecolor = :black,
            #markerstrokestyle = :dot,
            xlims = xlims,
            ylims = ( minimum(y) - ypad, maximum(y) + ypad ),
            ytickfontsize = 10,
            xtickfontsize = 10,
            xminorticks=minorticks,
            grid=:off,
            framestyle=:box,
            label=label1,
            legendfontsize=10,
            tickorientation=:in,
            legend=legend,
            size=(600,400)
          )
    
    if length(xtickvals) > 0
        plt = xticks!(xtickvals)
    end
    
    plt = yticks!(y1, n) 
    plt = xlabel!(xlabel, xguidefontsize=14)
    plt = vline!([1], linestyle=:dot, linewidth=1.5, linecolor=:black, label="")

    #Plot the second series if provided as m2, m2e1, m2e2
    if length(m2) > 0

        if length(m2) != length(m)
            error("Length of second data series must match first data series.")
        end
        
        if CI == true
            le2 = abs.(m2 .- m2e1)
            ue2 = abs.(m2e2 .- m2)
        else
            if length(e2) > 0
                println("Found values for m2e2 but setting m2e2 == e1 because CI != true")
            end
            le2 = m2e1
            ue2 = m2e1
        end

        err2 = Tuple.(zip(le2,ue2))
        
        # Don't scale axis on second plot!
        plt = scatter!(m2, y2, xaxis=scale_x, xerr=err2,
                 markersize = markersize,
                 markeralpha = 0.99,
                 markercolor = :firebrick4,
                 markerstrokewidth = 1,
                 markerstrokealpha = 0.99,
                 markerstrokecolor = :firebrick4,
                       label=label2,
                       legend_font_pointsize=9
                 )


        plt = yticks!(reverse(y), n)
        plt = xlabel!(xlabel, xguidefontsize=12)

        yc = y
        push!(yc, (y[end] + 1.0))
        pushfirst!(yc, 0)

        plt = hspan!([yc .- 0.5], alpha=0.3, color=:snow, label="")
        plt = hspan!([yc .+ 0.5], alpha=0.3, color=:grey80, label="")

    end

    display(plt)

end
