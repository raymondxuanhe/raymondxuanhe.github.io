#Julia code to create bar plots showing the CPU benchmarks
#Author: Raymond He
#Last updated: July 26 2026

using CairoMakie

# ── same palette as the bar charts ─────────────────────────────────────────
const BG      = "#0b2545"   # figure background
const PANEL   = "#0f2f5c"   # plot area
const C_M5    = "#5b9bd5"   # M5 Pro
const C_RYZEN = "#ed7d31"   # Ryzen 9 9950X
const FAINT   = ("#ffffff", 0.10)
const SUBTLE  = "#a9c6ea"

# ── data (seconds) ─────────────────────────────────────────────────────────
threads = 1:12

m5    = [112.5, 62.3, 44.1, 34.7, 28.9, 24.7, 21.9, 19.8, 17.9, 17.2, 15.9, 15.4]
ryzen = [140.0, 73.5, 50.6, 39.7, 33.7, 28.3, 25.5, 23.1, 21.3, 19.9, 18.6, 17.3]

fmt(v) = string(round(v; digits = 1))

# ── figure ─────────────────────────────────────────────────────────────────
fig = Figure(
    size = (1000, 560),
    backgroundcolor = BG,
    fonts = (; regular = "TeX Gyre Heros Makie",
               bold    = "TeX Gyre Heros Makie Bold"),
)

ax = Axis(fig[1, 1],
    title       = "Arellano (2008) — Time to Solve (CPU)",
    titlecolor  = :white,
    titlesize   = 26,
    titlefont   = :bold,
    titlegap    = 6,
    subtitle    = "Time to solve vs. thread count · lower is better",
    subtitlecolor = SUBTLE,
    subtitlesize  = 15,
    subtitlegap   = 12,

    backgroundcolor = PANEL,

    xlabel = "threads",
    ylabel = "seconds",
    xlabelcolor = SUBTLE, ylabelcolor = SUBTLE,
    xlabelsize  = 14,     ylabelsize  = 14,

    xticks = 1:12,
    yticks = 0:20:140,

    xticklabelcolor = :white, xticklabelsize = 13,
    yticklabelcolor = :white, yticklabelsize = 13,

    xgridcolor = FAINT, xgridwidth = 1,
    ygridcolor = FAINT, ygridwidth = 1,

    xtickcolor = ("#ffffff", 0.3), ytickcolor = ("#ffffff", 0.3),
    bottomspinecolor = ("#ffffff", 0.35),
    leftspinecolor   = ("#ffffff", 0.35),
    topspinevisible   = false,
    rightspinevisible = false,
)

# lines first, markers on top so the dots sit above the stroke
lines!(ax, threads, ryzen; color = C_RYZEN, linewidth = 3, label = "Ryzen 9 9950X")
lines!(ax, threads, m5;    color = C_M5,    linewidth = 3, label = "M5 Pro")

scatter!(ax, threads, ryzen; color = C_RYZEN, markersize = 11,
         strokecolor = PANEL, strokewidth = 1.5)
scatter!(ax, threads, m5;    color = C_M5,    markersize = 11,
         strokecolor = PANEL, strokewidth = 1.5)

# point labels: Ryzen above its curve, M5 below, so they never collide
# where the two series converge at high thread counts
text!(ax, threads, ryzen; text = fmt.(ryzen),
      align = (:center, :bottom), offset = (0, 11),
      color = :white, fontsize = 12, font = :bold)

text!(ax, threads, m5; text = fmt.(m5),
      align = (:center, :top), offset = (0, -11),
      color = :white, fontsize = 12, font = :bold)

xlims!(ax, 0.5, 12.5)
ylims!(ax, 0, 158)      # headroom for the label on the 140.0 point

axislegend(ax;
    position     = :rt,
    framevisible = false,
    labelcolor   = :white,
    labelsize    = 15,
    patchsize    = (26, 14),
    padding      = (8, 8, 8, 8),
)

fig

save(pwd()*"/files/benchmarks/cpu_benchmark.png", fig; px_per_unit = 4)

# ── same palette as the other charts ───────────────────────────────────────
const BG     = "#0b2545"   # figure background
const PANEL  = "#0f2f5c"   # plot area
const C_FP32 = "#5b9bd5"   # single precision
const C_FP64 = "#ed7d31"   # double precision
const FAINT  = ("#ffffff", 0.10)
const SUBTLE = "#a9c6ea"

# ── data (seconds) ─────────────────────────────────────────────────────────
# One row per card, two dodge slots for precision.  GPUs only, so the span is
# just ~14x — plain linear seconds, no log axis and no rebasing needed.
hardware = ["RTX 5090", "RTX 4080 Super"]

#       5090    5090    4080S   4080S
#       FP32    FP64    FP32    FP64
cat  = [1,      1,      2,      2     ]
secs = [0.12,   0.79,   0.32,   1.69  ]
grp  = [1,      2,      1,      2     ]
cols = [C_FP32, C_FP64, C_FP32, C_FP64]

fmt(v) = string(round(v; digits = 2), " s")

# ── figure ─────────────────────────────────────────────────────────────────
fig = Figure(
    size = (1000, 360),
    backgroundcolor = BG,
    fonts = (; regular = "TeX Gyre Heros Makie",
               bold    = "TeX Gyre Heros Makie Bold"),
)

ax = Axis(fig[1, 1],
    title       = "Arellano (2008) — Time to Solve (GPU)",
    titlecolor  = :white,
    titlesize   = 26,
    titlefont   = :bold,
    titlegap    = 6,
    subtitle    = "Single vs. double precision · lower is better",
    subtitlecolor = SUBTLE,
    subtitlesize  = 15,
    subtitlegap   = 12,

    backgroundcolor = PANEL,

    xlabel      = "seconds",
    xlabelcolor = SUBTLE,
    xlabelsize  = 14,

    yticks    = (1:length(hardware), hardware),
    yreversed = true,                       # row 1 at the top
    xticks    = 0:0.25:1.75,

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

barplot!(ax, cat, secs;
    direction  = :x,
    dodge      = grp,
    n_dodge    = 2,
    color      = cols,
    gap        = 0.32,
    dodge_gap  = 0.06,
    bar_labels = fmt.(secs),
    label_size   = 15,
    label_color  = :white,
    label_font   = :bold,
    label_offset = 6,
)

xlims!(ax, 0, 1.95)   # headroom for the end-of-bar labels

Legend(fig[2, 1],
    [PolyElement(color = c) for c in (C_FP32, C_FP64)],
    ["FP32", "FP64"];
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

save(pwd()*"/files/benchmarks/gpu_benchmark.png", fig; px_per_unit = 4)
