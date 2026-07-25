#Julia code to create bar plots showing the CPU benchmarks
#Author: Raymond He
#Last updated: July 25 2026


using CairoMakie

#################################################
####### Solving the Arellano (2008) model #######
#################################################
const BG       = "#0b2545"   # figure background
const PANEL    = "#0f2f5c"   # plot area
const C_FP64   = "#3fa89c"   # GPU, double precision
const C_MULTI  = "#5b9bd5"   # CPU, multi-threaded
const C_SINGLE = "#ed7d31"   # CPU, single-threaded
const FAINT    = ("#ffffff", 0.10)
const SUBTLE   = "#a9c6ea"

const LOGSCALE = false        # flip to false for a linear relative-time axis

# ── data ───────────────────────────────────────────────────────────────────
# FP32 is the yardstick, not a bar: bars are plotted as multiples of it,
# while the white end labels stay in raw seconds.
const BASELINE = 0.12        # RTX 5090, FP32, seconds

# ordered fastest → slowest so the winner sits at the top
hardware = ["RTX 5090", "Ryzen 9 9950X", "M4 Max"]

cat  = [1,      2,       3,       2,        3       ]   # row index
secs = [0.79,   17.2,    24.7,    140.0,    152.1   ]   # raw seconds
grp  = [1,      1,       1,       2,        2       ]   # dodge group
cols = [C_FP64, C_MULTI, C_MULTI, C_SINGLE, C_SINGLE]

rel = secs ./ BASELINE       # what actually gets drawn

fmt(v) = v < 10 ? string(round(v; digits = 2), " s") : string(round(v; digits = 1), " s")

# ── scale-dependent settings ───────────────────────────────────────────────
# On a log axis 1x is the multiplicative origin, so bars grow from the
# baseline itself and the left spine *is* the FP32 reference.
scale_kwargs = LOGSCALE ?
    (; xscale = log10,
       xticks = ([1, 10, 100, 1000], ["1×", "10×", "100×", "1,000×"])) :
    (; xticks = (0:200:1400, string.(0:200:1400, "×")))

fillto = LOGSCALE ? 1.0 : 0.0
xrange = LOGSCALE ? (1.0, 3000) : (0, 1500)

# ── figure ─────────────────────────────────────────────────────────────────
fig = Figure(
    size = (1000, 430),
    backgroundcolor = BG,
    fonts = (; regular = "TeX Gyre Heros Makie",
               bold    = "TeX Gyre Heros Makie Bold"),
)

ax = Axis(fig[1, 1];
    title       = "Arellano (2008) — Time to Solve",
    titlecolor  = :white,
    titlesize   = 26,
    titlefont   = :bold,
    titlegap    = 6,
    subtitle    = "Relative to RTX 5090 FP32 (0.12 s = 1×) · lower is better",
    subtitlecolor = SUBTLE,
    subtitlesize  = 15,
    subtitlegap   = 12,

    backgroundcolor = PANEL,

    xlabel      = LOGSCALE ? "slowdown vs. FP32 baseline (log scale)" :
                             "slowdown vs. FP32 baseline",
    xlabelcolor = SUBTLE,
    xlabelsize  = 14,

    yticks    = (1:length(hardware), hardware),
    yreversed = true,                       # row 1 at the top

    xticklabelcolor = :white, xticklabelsize = 13,
    yticklabelcolor = :white, yticklabelsize = 16,

    xgridcolor    = FAINT, xgridwidth = 1,
    ygridvisible  = false,
    yticksvisible = false,
    xtickcolor    = ("#ffffff", 0.3),

    bottomspinecolor = ("#ffffff", 0.35),
    leftspinecolor   = ("#ffffff", 0.55),   # this spine is the baseline
    topspinevisible   = false,
    rightspinevisible = false,

    scale_kwargs...
)

barplot!(ax, cat, rel;
    direction  = :x,
    dodge      = grp,
    n_dodge    = 2,
    color      = cols,
    fillto     = fillto,
    gap        = 0.30,
    dodge_gap  = 0.06,
    bar_labels = fmt.(secs),               # white labels stay in seconds
    label_size   = 15,
    label_color  = :white,
    label_font   = :bold,
    label_offset = 6,
)

xlims!(ax, xrange...)   # headroom for the end-of-bar labels

Legend(fig[2, 1],
    [PolyElement(color = c) for c in (C_FP64, C_MULTI, C_SINGLE)],
    ["GPU — FP64", "CPU — multi-threaded", "CPU — single-threaded"];
    orientation   = :horizontal,
    framevisible  = false,
    labelcolor    = :white,
    labelsize     = 14,
    tellheight    = true,
    tellwidth     = false,
    padding       = (0, 0, 0, 0),
)
rowgap!(fig.layout, 8)

fig

save(pwd()*"/files/benchmarks/arellano_2008.png", fig; px_per_unit = 4)



#####################################################
####### Cleaning & Processing A Large Dataset #######
#####################################################
const BG       = "#0b2545"   # figure background
const PANEL    = "#0f2f5c"   # plot area
const C_MULTI  = "#5b9bd5"   # multi-threaded
const C_SINGLE = "#ed7d31"   # single-threaded
const FAINT    = ("#ffffff", 0.10)
const SUBTLE   = "#a9c6ea"

# ── data ───────────────────────────────────────────────────────────────────
# span is only ~5.8x here, so a linear axis is fine — no log, no rebasing
hardware = ["Ryzen 9 9950X", "M4 Max"]

cat  = [1,       2,       1,        2       ]   # row index
mins = [59.9,    79.0,    152.5,    346.4   ]   # minutes
grp  = [1,       1,       2,        2       ]   # dodge group
cols = [C_MULTI, C_MULTI, C_SINGLE, C_SINGLE]

fmt(v) = string(round(v; digits = 1), " min")

# ── figure ─────────────────────────────────────────────────────────────────
fig = Figure(
    size = (1000, 360),
    backgroundcolor = BG,
    fonts = (; regular = "TeX Gyre Heros Makie",
               bold    = "TeX Gyre Heros Makie Bold"),
)

ax = Axis(fig[1, 1],
    title       = "Cleaning & Processing A Large Dataset - Wall Clock",
    titlecolor  = :white,
    titlesize   = 26,
    titlefont   = :bold,
    titlegap    = 6,
    subtitle    = "Full clean-and-process pass · lower is better",
    subtitlecolor = SUBTLE,
    subtitlesize  = 15,
    subtitlegap   = 12,

    backgroundcolor = PANEL,

    xlabel      = "minutes",
    xlabelcolor = SUBTLE,
    xlabelsize  = 14,

    yticks    = (1:length(hardware), hardware),
    yreversed = true,                       # row 1 at the top
    xticks    = 0:50:350,

    xticklabelcolor = :white, xticklabelsize = 13,
    yticklabelcolor = :white, yticklabelsize = 16,

    xgridcolor    = FAINT, xgridwidth = 1,
    ygridvisible  = false,
    yticksvisible = false,
    xtickcolor    = ("#ffffff", 0.3),

    bottomspinecolor = ("#ffffff", 0.35),
    leftspinecolor   = ("#ffffff", 0.35),
    topspinevisible   = false,
    rightspinevisible = false,
)

barplot!(ax, cat, mins;
    direction  = :x,
    dodge      = grp,
    n_dodge    = 2,
    color      = cols,
    gap        = 0.32,
    dodge_gap  = 0.06,
    bar_labels = fmt.(mins),
    label_size   = 15,
    label_color  = :white,
    label_font   = :bold,
    label_offset = 6,
)

xlims!(ax, 0, 395)   # headroom for the end-of-bar labels

Legend(fig[2, 1],
    [PolyElement(color = c) for c in (C_MULTI, C_SINGLE)],
    ["Multi-threaded", "Single-threaded"];
    orientation   = :horizontal,
    framevisible  = false,
    labelcolor    = :white,
    labelsize     = 14,
    tellheight    = true,
    tellwidth     = false,
    padding       = (0, 0, 0, 0),
)
rowgap!(fig.layout, 8)

save(pwd()*"/files/benchmarks/lightcast.png", fig; px_per_unit = 4)
fig
