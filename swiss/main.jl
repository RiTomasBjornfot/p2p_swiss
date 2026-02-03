#=
    Samma lösning som i Schweiz för TCC, ICC, HNA/LNA:
    TCC -> FL1-A, FL2-A
    ICC -> FL1-A, FL2-A
    HNA -> FL1-A
=#

using DelimitedFiles
using TOML
using PolygonOps, StaticArrays
using GLMakie
import Pipe.@pipe as @p

#   ==================
#       DATA IMPORT
#   ==================

"""
    read_data(conf, fname) 

    reads the channels defines in configuration a matrix from csv file
    conf:   The configureation as a dictionary
    fname:  The name of the file
    Returns the log10(x+1) of all values
"""
function read_data(conf, fname) 
    x, header = readdlm(conf["dir"]*fname, ',', header=true)
    X = zeros(size(x, 1), length(conf["channels"]))
    for (i, h) in enumerate(conf["channels"])
        j = findfirst(x -> x == h, vec(header))
        X[:, i] = x[:, j]
    end
    @. log10(X + 1)
end

"""
    index = get_channel_index(conf, channel)
    Finds the index of a channel 
    conf:       The configureation as a dictionary
    channel:    The channel to find
    Returns the index
"""
get_channel_index(conf, channel) = findfirst(c -> c == channel, conf["channels"])

#   ==============
#       GATES
#   ==============

"""
    y = tgate(x, op, gate)
    
    Counts events for a 1D gate
    x:          the data (FL1-A, FL2-A or SSC-A)
    op:         the operator (> or <)
    gate:       the gate (threshold)      
"""
tgate(x, op, gate) = sum(op.(x, gate))

"""
    poly = gate2poly(gate)

    Converts the gate to a polygon according to PolygonOps
    X:      The data as a nX2 matrix
    lim:   The gate 
    Returns the values as nX2 matrix
"""

gate2poly(gate) = @p map(i -> append!(gate[i], gate[i][1]), 1:2) |> SVector.(_...)

"""
    Y = pgate(X, gate)

    Finds the values inside a polygon.
    X:      The channel values as a nX2 matrix
    gate:   The polygon gate
    Returns a nX2 matrix
"""
function pgate(X, gate)
    poly, j = gate2poly(gate), []
    for i in axes(X, 1)
        p = SVector(X[i, :]...)
        if inpolygon(p, poly) == 1
            push!(j, i)
        end
    end
    X[j, :]
end

"""
    bin(x, lim, sz)

    Put the data in bins. 
    z:      The data vector
    lim:    The upper and lower bounds
    sz:     The number of bins
    
    Return: the bin-number each element in the vector belongs to.
    NOTE: All elements in z must be inside the boundary!
"""
function bin(x, lim, sz)
    x = x .- lim[1]
    x = x ./ (lim[2] - lim[1])
    x = x*sz .+ 0.5
    round.(Int, x)
end

"""
    fingerprint(conf, X, channel)

    Makes a 1D fingerprint
    conf:       The configuration
    X:          All data as a matrix where each column is a channel (e.g. FL1-A, FL2-A or SSC-A)
    channel:    The channel to use
"""
function fingerprint(conf, X, channel)
    i = get_channel_index(conf, channel)
    x = X[:, i] 
    lim = conf["fp_bound"][i]
    sz = conf["fp_sz"][i]
    x = bin(x, lim, sz)
    map(n -> sum(x .== n), 1:sz)
end

"""
    Z = fingerprint2d(conf, X, channels)

    Makes a 2D fingerprint
    conf:       The configuration
    X:          All data as a matrix where each column is a channel (e.g. FL1-A, FL2-A or SSC-A)
    channels:   The channels to use
"""
function fingerprint2d(conf, X, channels)
    # get the acctual data rows
    i, j = map(k -> get_channel_index(conf["import"], channels[k]), 1:2)
    Z = X[:, [i, j]]
    # get the boundries and sizes
    lim = @p conf["fingerprint"]["boundary"] |> map(k -> _[k], [i, j])
    sz = @p conf["fingerprint"]["size"] |> map(k -> _[k], [i, j])

    # removes all events outside fp_bound
    Z = @p findall(z -> lim[1][1] < z < lim[1][2], Z[:, 1]) |> Z[_, :]
    Z = @p findall(z -> lim[2][1] < z < lim[2][2], Z[:, 2]) |> Z[_, :]
    
    # binning
    x, y = map(i -> bin(Z[:, i], lim[i], sz[i]), 1:2)
    mat = zeros(sz...)
    for i in axes(x, 1)
        mat[x[i], y[i]] += 1
    end
    mat'
end

"""
    plot_fingerprint2d(mat, channels)

    A nice plot of the 2D fingerprint
    mat:        The fingerprint (from fingerprint2d)
    channels:   The channels used
    Returns a figure
"""
function plot_fingerprint2d(mat, channels)
    fig = Figure()
    ax = Makie.Axis(fig[1, 1])
    hm = heatmap!(ax, mat', interpolate=false, ; colormap=:viridis)
    Colorbar(fig[1, 2], hm, label="Events")
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.xlabel = channels[1]
    ax.ylabel = channels[2]
    fig
end

"""
    plot_gate(X, gate)

    A nice plot of the gate and scatter data
    X:      The events
    gate:   The gate as a polygon
    Returns a figure
"""
function plot_gate(X, gate)
    fig = Figure()
    ax = Makie.Axis(fig[1, 1])
    scatter!(ax, X, markersize=3)
    x = append!(gate[1], gate[1][1])
    y = append!(gate[2], gate[2][1])
    lines!(ax, x, y, color=:red)
    fig
end

#   ============
#       MAIN    
#   ============

conf = TOML.parsefile("config.toml")

X = read_data(conf["import"], "Sample_A01.csv")
println("Total number of events: ", size(X, 1))

# TCC
X_tcc = pgate(X, conf["gate"]["tcc"])
println("Total Cell Count (TCC): ", size(X_tcc, 1))

# ICC
X_icc = pgate(X, conf["gate"]["icc"])
println("Intact Cell Count (ICC): ", size(X_icc, 1))

# HNA
no_hna = tgate(X_icc[:, 1], >, conf["gate"]["hna"])
println("Number of HNA events: ", no_hna)

fig_tcc =plot_gate(X, conf["gate"]["tcc"])

fp = fingerprint2d(conf, X_tcc, ["FL1-A", "FL2-A"])
fig_fp = plot_fingerprint2d(fp, ["FL1-A", "FL2-A"])
