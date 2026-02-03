using DelimitedFiles
using Images
using BSplineKit
using Peaks
using Statistics
using PolygonOps
using StaticArrays
using GLMakie
import Pipe.@pipe as @p
import Contour as cont


# the 10th log is used in the calculation
tolog(x) = @. log10(x + 1)

# get the channel vectors from the csv 
function read_csv(fp, channels)
    data = @p split(fp, ".")[1]*".csv" |> readdlm(_, ',')
    X = []
    for ch ∈ channels
        x = @p findfirst(x -> x == ch, data[1, :]) |> data[2:end, _] |> tolog
        push!(X, x)
    end
    X
end

# Bins the x vector in into n bins 
function binning(x, n)
    z = x .- minimum(x)    
    z ./= maximum(z)
    z .*= (n - 1)
    round.(Int, z) .+ 1
end

# makes a 2D histogram of x and y vectors
function bin_img(x, y, n)
    xb, yb = binning(x, n), binning(y, n)
    X = zeros(Int, (n, n))
    for i ∈ 1:length(xb)
        X[yb[i], xb[i]] += 1 
    end
    X
end

# Calculates the bin center values and bin edges
# x is the vector to be binned
function bin_limits(x, n)
    bin_edges = range(minimum(x), maximum(x), n + 1)
    bin_cc = Float64[]
    for i ∈ 1:(length(bin_edges)-1)
        cc = (bin_edges[i] + bin_edges[i+1]) / 2
        push!(bin_cc, cc)
    end
    bin_cc, bin_edges
end

# Closes a polygon if not closed.
close_poly(x) = isequal(x[1, :], x[end, :]) ? x : vcat(x, [x[1, 1] x[1, 2]])

# Calculates the area of a contour
carea(x) = SVector.(x[:, 1], x[:, 2]) |> PolygonOps.area |> abs

# Finds contours in a binary image
function cnts(z, n)
    z = Float64.(z)
    cnt, sz = Matrix{Float64}[], size(z)
    for cl ∈ cont.levels(cont.contours(1:sz[1], 1:sz[2], z, n))
        println(length(cont.lines(cl)))
        for line ∈ cont.lines(cl)
            α, β  = cont.coordinates(line)
            push!(cnt, [α β])
        end
    end
    close_poly.(cnt)
end

# kollar om några ligger innanför någon annan av konturerna
function outer_contours(c)
    # tar ut första punkten från alla c
    b = zeros(Bool, length(c))
    for i ∈ axes(c, 1)
        k = size(c[i], 1) ÷ 2
        p = SVector(c[i][k, :]...)
        for j ∈ axes(c, 1)
            (j == i) && (continue)
            poly = SVector.(c[j][:, 1], c[j][:, 2])
            if inpolygon(p, poly) == 1
                b[i] = true
                break
            end
        end
    end
    b_outside = findall(b₀ -> b₀ == false, b)
    c[b_outside]
end

# Räknar ut en kurva på topparna av histogrammet
function bin_filter(sx, y)
    n = length(y)
    xf = interpolate(1:n, y, BSplineOrder(2))
    imfilter(xf.(sx), Kernel.gaussian((11,)))
end

#=
#  "ICC Gaten" som ligger mellan intakta och inte intakta celler defineras som:
#   * Ett minima mellan dom två största maximan på FL3-A histogrammet.
#       - Om det inte finns två maximan ges ingen gate (error?).
#   * Om det finns fler minima mellan maximam, är "gaten" medianen av miniman
=#
function calc_icc_gate(sx, sy)
    mx = findmaxima(sy, 10)
    # If less than two maximas are found no gate is given
    (length(mx.indices) < 2) && (return -1)
    # if more than two maximas if found we take the two highest ones
    i = sortperm(mx.heights, rev=true)[1:2]
    pks_index = round.(Int, mx.indices[i])
    mn = findminima(sy, 10)
    # finds the minimas between the maximas
    i = findall(mn -> minimum(pks_index) < mn < maximum(pks_index), mn.indices)
    # if no minima is found, the gate is the mean of the maximas
    # if minimas are found inside the maximas the gate is the median of indices
    k = length(i) == 0 ? sum(pks_index)/2 : median(mn.indices[i])
    gate_index = round(Int, k)
    gate_value = sx[gate_index]
    return round(Int, gate_value)
end

rnd(x) = round(x, digits=2)


# plottar histogramlösningen för icc
function icc_hist_plot(X2, sx, sy, gate)
    fig = Figure(size=(1000, 500))
    ax1 = Makie.Axis(fig[1, 1], ylabel=channels[2], xlabel=channels[1], title=fp)
    ax1.xticksvisible = false
    ax1.xticklabelsvisible = false
    ax1.yticksvisible = false
    ax1.yticklabelsvisible = false
    image!(ax1, X2', interpolate=false)
    hlines!(ax1, gate, color=:red)

    ax2 = Makie.Axis(fig[1, 2], xlabel="Events")
    barplot!(ax2, y3, direction=:x)
    lines!(ax2, sy, sx)
    hlines!(ax2, gate, color=:red)
    fig
end


# ========
#   MAIN
# ========
fp = "testdata/Sample_A05.fcs"
n = 50

# hittar på gates
fl1_gate = 2.9
ssc_gate = 0.2

channels = ["FL1-A", "FL3-A", "SSC-A"]

# läser in csv filen 
if !isfile(split(fp, ".")[1]*".csv")
     #läser in fcs filen och konverterar till csv
    `.venv/bin/python3 fcs_convert.py $fp` |> run
end
Z_org = read_csv(fp, channels)
sum_of_events = length(Z_org[1])

# gör data till bilder
X1 = bin_img(Z_org[1], Z_org[3], n) # för TCC och HNA/LNA
X2 = bin_img(Z_org[1], Z_org[2], n) # för ICC

# Tar ut bin center och kanter
fl1_cc, fl1_eg = bin_limits(Z_org[1], n)
ssc_cc, ssc_eg = bin_limits(Z_org[2], n)
fl3_cc, fl3_eg = bin_limits(Z_org[3], n)

# --- Steg 1: Noise Reduction ---
# NOTE: Alla värden som ligger i samma bin som gaten tas bort!
fl1_bin_gate = findfirst(x -> x > fl1_gate, fl1_eg) - 1
ssc_bin_gate = findfirst(x -> x > ssc_gate, ssc_eg) - 1

tcc_fig = Figure(size=(1000, 500))
ax = Makie.Axis(tcc_fig[1, 1], xlabel=channels[1], ylabel=channels[3])
ax.title = "All events"
image!(ax, X1', interpolate=false)
vlines!(ax, fl1_bin_gate, color=:red)
hlines!(ax, ssc_bin_gate, color=:red)

# Här tas bruset bort!
if fl1_bin_gate > 0
    X1[:, 1:fl1_bin_gate] .= 0
end
if ssc_bin_gate > 0
    X1[1:ssc_bin_gate, :] .= 0
end

ax2 = Makie.Axis(tcc_fig[1, 2], xlabel=channels[1], ylabel=channels[3])
ax2.title = "Noise removed"
image!(ax2, X1', interpolate=false)
vlines!(ax2, fl1_bin_gate, color=:red)
hlines!(ax2, ssc_bin_gate, color=:red)
linkaxes!(ax, ax2)
save("out/tcc_fig.png", tcc_fig)

# --- Räknar ut TCC ---
tcc = sum(X1)

# printar TCC
println('-'^20)
println("No events: ", length(Z_org[1]))
println("TCC: ", tcc)


# --- Räknar ut ICC  med minimat på histogrammet ---
y3 = vec(sum(X2, dims=2))
sx = range(1, n, 10*n)
sy = bin_filter(sx, y3)
icc_gate = calc_icc_gate(sx, sy)
icc_gate_value = fl3_cc[icc_gate]

icc = sum(X2[:, 1:icc_gate])
icc_precent = 100*icc/sum(X2)

println('-'^20)
println("ICC:\t", icc)
println("% ICC:\t", rnd(icc_precent))
println("Gate: ", rnd(icc_gate_value))
println('-'^20)

# plottar ICC "hist" beräkningen
icc_hist_fig = icc_hist_plot(sqrt.(X2), sx, sy, icc_gate)
save("out/icc_hist_fig.png", icc_hist_fig)

# Räknar ut ICC med clustering

# --- Räknar ut ICC med conturer ---
# tar ut konturer på 50 nivåer
c = cnts(X2, 50)
# behåller bara dom yttre konturerna 
c = outer_contours(c)
# behåller bara dom två största klustren
c = @p carea.(c) |> sortperm(_, rev=true)[1:2] |> c[_]

icc_cluster_fig = Figure()
ax3 = Makie.Axis(icc_cluster_fig[1, 1])
image!(ax3, sqrt.(X2)', interpolate=false)
for c₀ ∈ c
    lines!(ax3, c₀[:, 2], c₀[:, 1], color=:red)
end
save("out/icc_cluster_fig.png", icc_cluster_fig)
