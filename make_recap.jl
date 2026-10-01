using ArgParse, CairoMakie

# Assemble the per-simulation figures written by exec_analysis.jl into recap grids:
# one PNG per figure kind (polar, cluster, rdf, situa), rows are potential × packing
# fraction, columns are offcenter.

# ── Command-line arguments ─────────────────────────────────────────────────────
let s = ArgParseSettings(description = "Recap grids of the exec_analysis.jl figures.")
    @add_arg_table! s begin
        "imgs_dir"
            help     = "folder with the exec_analysis.jl figures"
            required = true
        "out_dir"
            help     = "folder for the recap PNGs"
            required = true
        "--label"
            help     = "title suffix, e.g. \"run02, range\""
            required = true
        "--prefix"
            help     = "file name prefix (default: last word of label)"
            default  = nothing
    end
    global args = parse_args(s)
end

imgs_dir = args["imgs_dir"]
out_dir  = args["out_dir"]
label    = args["label"]
prefix   = something(args["prefix"], strip(last(split(label, ","))))

const POTENTIALS = [("lj_nondiv", "LJ nondiv"),
                    ("scaled_attractive", "attractive"),
                    ("scaled_repulsive", "repulsive")]
const PFS        = [("1e-01", "0.1"), ("2e-01", "0.2")]
const OFFCENTERS = [("-5e-01", "-0.5"), ("0e+00", "0"), ("5e-01", "+0.5")]
const KINDS      = ["polar", "cluster", "rdf", "situa"]

# Layout in pixels - matches the run01 recaps
const WIDTH    = 2586
const LABEL_W  = 253    # left column holding the row labels
const RIGHT_M  = 23
const HEADER_H = 184    # suptitle + column headers
const BOTTOM_M = 34
const FONT     = "DejaVu Sans Bold"

img_path(kind, p, f, o) = joinpath(imgs_dir, "$(kind)_$(p)_pf=$(f)_offcenter=$(o).png")

function make_recap(kind, out_path)
    rows  = [(p, pl, f, fl) for (p, pl) in POTENTIALS for (f, fl) in PFS]
    paths = [img_path(kind, p, f, o) for (p, _, f, _) in rows for (o, _) in OFFCENTERS]
    first_img = findfirst(isfile, paths)
    first_img === nothing && (@warn "no $kind figures in $imgs_dir"; return)
    h0, w0 = size(CairoMakie.load(paths[first_img]))

    cell_w = (WIDTH - LABEL_W - RIGHT_M) / length(OFFCENTERS)
    cell_h = cell_w * h0 / w0
    height = round(Int, HEADER_H + length(rows) * cell_h + BOTTOM_M)
    fy(px) = height - px                  # pixel from top -> Makie y (from bottom)

    fig = Figure(size = (WIDTH, height), figure_padding = 0, backgroundcolor = :white)
    text!(fig.scene, WIDTH / 2, fy(65), text = "$(uppercase(kind))  —  recap ($label)",
        font = FONT, fontsize = 52.8, align = (:center, :center))
    for (j, (_, ol)) in enumerate(OFFCENTERS)
        text!(fig.scene, LABEL_W + (j - 0.5) * cell_w, fy(151), text = "offcenter = $ol",
            font = FONT, fontsize = 29.2, align = (:center, :center))
    end

    for (i, (p, pl, f, fl)) in enumerate(rows)
        top = HEADER_H + (i - 1) * cell_h
        text!(fig.scene, 69, fy(top + cell_h / 2 - 5), text = "$pl\npf=$fl",
            font = FONT, fontsize = 26.7, align = (:left, :center), justification = :left)
        for (j, (o, _)) in enumerate(OFFCENTERS)
            left = LABEL_W + (j - 1) * cell_w
            ax = Axis(fig, bbox = BBox(left, left + cell_w, fy(top + cell_h), fy(top)))
            hidedecorations!(ax); hidespines!(ax)
            path = img_path(kind, p, f, o)
            if isfile(path)
                img = CairoMakie.load(path)
                image!(ax, (0, size(img, 2)), (0, size(img, 1)), rotr90(img))
                limits!(ax, 0, size(img, 2), 0, size(img, 1))
            else
                text!(ax, 0.5, 0.5, text = "missing", color = :gray, align = (:center, :center), space = :relative)
            end
        end
    end

    CairoMakie.save(out_path, fig, px_per_unit = 1)
    @info "wrote $out_path"
end

mkpath(out_dir)
for kind in KINDS
    make_recap(kind, joinpath(out_dir, "recap_$(prefix)_$(kind).png"))
end
