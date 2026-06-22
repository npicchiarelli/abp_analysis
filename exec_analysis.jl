using ArgParse, CairoMakie, Clustering, CSV, DataFrames, Distributions, JLD2, StatsBase
include("src/boundary_conditions.jl")
include("src/orderparameters.jl")
include("src/plot_animation.jl")
include("src/radialdistribution.jl")

# ── Command-line arguments ─────────────────────────────────────────────────────
let s = ArgParseSettings(description = "Cluster and polarization analysis for ABP simulations.")
    @add_arg_table! s begin
        "path"
            help     = "path to the run directory"
            required = true
        "--nbins"
            help     = "number of bins for the radial distribution function"
            arg_type = Int
            default  = 1000
        "--mindist-cluster"
            help     = "DBSCAN neighbourhood radius [μm]"
            arg_type = Float64
            default  = 5.0
        "--tlim"
            help     = "upper time limit for plots [s] (default: auto)"
            arg_type = Float64
            default  = nothing
    end
    global args = parse_args(s)
end

path            = args["path"]
nbins           = args["nbins"]
mindist_cluster = args["mindist-cluster"]
tlim            = args["tlim"]

# Parse a single key from a settings.jl file without executing it
function read_setting(settings_file, key)
    for line in eachline(settings_file)
        parts = split(line, '=', limit=2)
        if length(parts) == 2 && strip(parts[1]) == key
            return strip(first(split(parts[2], '#')))
        end
    end
    return nothing
end

# ── Result collectors (one entry per simulation directory) ─────────────────────
maximum_polarization = Float64[]
maxpoltime           = Float64[]
maxpolloc            = Float64[]
maxpolloc_time       = Float64[]
firstpeak_height     = Float64[]
fpk_pos              = Float64[]
oc                   = Float64[]

# ── Main loop over simulation directories ─────────────────────────────────────
for sim_dir in readdir(path, join = true)
    isdir(sim_dir) || continue
    settings_file = joinpath(sim_dir, "settings.jl")
    isfile(settings_file) || continue

    # Load simulation parameters from settings file
    L         = parse(Float64, read_setting(settings_file, "L"))
    R         = parse(Float64, read_setting(settings_file, "R"))
    Np        = parse(Int,     read_setting(settings_file, "Np"))
    offcenter = parse(Float64, read_setting(settings_file, "offcenter"))
    δt        = parse(Float64, read_setting(settings_file, "δt"))
    v         = parse(Float64, read_setting(settings_file, "v"))
    push!(oc, offcenter)

    df_list        = DataFrame[]
    df_meanor_list = DataFrame[]
    df_t_list      = DataFrame[]
    df_nclust_list = DataFrame[]

    # ── Inner loop over independent runs ──────────────────────────────────────
    for run_dir in readdir(sim_dir, join = true)
        isdir(run_dir) || continue
        run_name  = basename(run_dir)
        data_file = joinpath(run_dir, "simulation_$(run_name).txt")
        isfile(data_file) || continue

        df = CSV.read(data_file, DataFrame, ntasks = 16)
        push!(df_list, df)

        # Assign DBSCAN cluster labels, then correct for periodic boundaries
        df_cl = transform(
            groupby(df, :Time),
            [:xpos, :ypos] =>
                ((x,y) -> assignments(dbscan([x y]', mindist_cluster, min_cluster_size = 2))) =>
                :dbscan)
        periodic_clustering!(df_cl, L, mindist_cluster)

        # Cluster size statistics per timestep
        df_nclust = combine(
            groupby(df_cl[df_cl.dbscan .!== 0, :], [:Time, :dbscan]),
            nrow => :partcount)
        df_nclust = combine(
            groupby(df_nclust, :Time),
            :partcount => maximum => :max_cl_size,
            nrow => :nclust)

        # Local polarization: mean over clusters
        df_or = combine(
            groupby(df_cl, [:Time, :dbscan]),
            :orientation => mean_polarization => :polar,
            nrow => :partcount)
        df_meanor = combine(
            groupby(df_or[df_or.dbscan .!== nothing, :], [:Time]),
            :polar => mean => :mean_polar)

        # Global polarization: mean over all particles
        df_t = combine(groupby(df, :Time), :orientation => mean_polarization => :polar)

        push!(df_meanor_list, df_meanor)
        push!(df_t_list, df_t)
        push!(df_nclust_list, df_nclust)
    end

    # Snapshot of the last timestep of the first run
    plot_one_timestep(
        df_list[1], R, L, Int(size(df_list[1], 1) / Np),
        savepath = joinpath(path, "imgs", "situa_$(basename(sim_dir)).png"),
        display_plot = false)

    # Average observables across runs
    df_meannclust = combine(
        groupby(vcat(df_nclust_list...), :Time),
        :max_cl_size => mean => :max_cl_size,
        :nclust      => mean => :nclust,
        :max_cl_size => std  => :max_cl_size_std,
        :nclust      => std  => :nclust_std)
    df_meanor = combine(
        groupby(vcat(df_meanor_list...), :Time),
        :mean_polar => mean => :mean_polar,
        :mean_polar => std  => :std_polar)
    df_t = combine(
        groupby(vcat(df_t_list...), :Time),
        :polar => mean => :mean_polar,
        :polar => std  => :std_polar)

    # ── Transient time: first time polarization exceeds 95% of its maximum ────
    maxlocpol    = maximum(df_meanor.mean_polar)
    t_maxlocpol  = df_meanor.Time[findfirst(x -> x > 0.95*maxlocpol, df_meanor.mean_polar)] * δt
    println("local polar $maxlocpol")
    println("at time $t_maxlocpol")
    push!(maxpolloc, maxlocpol)
    push!(maxpolloc_time, t_maxlocpol)

    maxpol   = maximum(df_t.mean_polar)
    t_maxpol = df_t.Time[findfirst(x -> x > 0.95*maxpol, df_t.mean_polar)] * δt
    println("maximum polarization $maxpol")
    println("at time $t_maxpol")
    push!(maximum_polarization, maxpol)
    push!(maxpoltime, t_maxpol)

    # ── Plot: local and global polarization vs time ───────────────────────────
    p2 = Figure(fontsize = 20)
    ax = Axis(p2[1,1],
        limits = ((0, tlim), nothing),
        xlabel = "Time [s]", ylabel = "Polarization",
        xgridvisible = false, ygridvisible = false)
    lines!(ax, df_meanor.Time*δt, df_meanor.mean_polar, linewidth = 1., color = :blue,  label = "local")
    band!(ax,  df_meanor.Time*δt,
        df_meanor.mean_polar .- df_meanor.std_polar,
        df_meanor.mean_polar .+ df_meanor.std_polar,
        color = :blue, alpha = 0.3)
    lines!(ax, df_t.Time*δt, df_t.mean_polar, linewidth = 1., color = :green, label = "global")
    band!(ax,  df_t.Time*δt,
        df_t.mean_polar .- df_t.std_polar,
        df_t.mean_polar .+ df_t.std_polar,
        color = :green, alpha = 0.3)
    axislegend(ax)
    CairoMakie.save(joinpath(path, "imgs", "polar_$(basename(sim_dir)).png"), p2)

    # ── Plot: cluster size and count vs time ──────────────────────────────────
    p3  = Figure(fontsize = 20)
    ax1 = Axis(p3[1,1],
        limits = ((0, tlim), nothing),
        xlabel = "Time [s]", ylabel = "Max cluster size",
        xgridvisible = false, ygridvisible = false)
    ax2 = Axis(p3[2,1],
        limits = ((0, tlim), nothing),
        xlabel = "Time [s]", ylabel = "Number of clusters",
        xgridvisible = false, ygridvisible = false)
    linkxaxes!(ax1, ax2)
    lines!(ax1, df_meannclust.Time*δt, df_meannclust.max_cl_size, linewidth = 2., color = :black)
    band!(ax1,  df_meannclust.Time*δt,
        df_meannclust.max_cl_size .- df_meannclust.max_cl_size_std,
        df_meannclust.max_cl_size .+ df_meannclust.max_cl_size_std,
        color = :black, alpha = 0.3)
    lines!(ax2, df_meannclust.Time*δt, df_meannclust.nclust, linewidth = 2., color = :orange)
    band!(ax2,  df_meannclust.Time*δt,
        df_meannclust.nclust .- df_meannclust.nclust_std,
        df_meannclust.nclust .+ df_meannclust.nclust_std,
        color = :orange, alpha = 0.3)
    CairoMakie.save(joinpath(path, "imgs", "cluster_$(basename(sim_dir)).png"), p3)

    # ── Plot: radial distribution function at the last timestep ───────────────
    rdf = combine(
        groupby(df_list[1], :Time),
        [:xpos, :ypos] =>
            ((x,y) -> (radialdistributionfunction(copy(x), copy(y), R, L, nbins),)) =>
            :RadialDistributionFunction)
    transform!(rdf, :RadialDistributionFunction => ByRow(x -> x[1]) => :RadialDistributionFunction)
    rs, bins = radialbinssquare(L, nbins)

    p4 = Figure(fontsize = 20)
    ax = Axis(p4[1,1],
        xlabel = "Distance [μm]", ylabel = "Radial Distribution Function",
        xgridvisible = false, ygridvisible = false)
    lines!(ax, rs, rdf.RadialDistributionFunction[end], linewidth = 2.)

    firstpeak     = maximum(rdf.RadialDistributionFunction[end])
    firstpeak_pos = rs[argmax(rdf.RadialDistributionFunction[end])]
    println("first peak $firstpeak")
    println("fp position $firstpeak_pos")
    push!(firstpeak_height, firstpeak)
    push!(fpk_pos, firstpeak_pos)
    CairoMakie.save(joinpath(path, "imgs", "rdf_$(basename(sim_dir)).png"), p4)
end

# ── Summary: transient time vs offcenter across all simulations ────────────────
fig = Figure(fontsize = 20)
ax  = Axis(fig[1,1],
    xgridvisible = false, ygridvisible = false,
    xlabel = "offcenter", ylabel = "Transient Time [s]")
CairoMakie.scatter!(ax, oc, maxpoltime, markersize = 20, color = :black)
CairoMakie.save(joinpath(path, "imgs", "transient_time.png"), fig)
