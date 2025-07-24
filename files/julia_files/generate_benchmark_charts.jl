#Julia code to create bar plots showing the CPU benchmarks
#Author: Raymond He
#Last updated: June 5 2025

using Plots

#################################################
####### Solving the Arellano (2008) model #######
#################################################

# CPU performance data
cpus = ["Ryzen 9 9950X", "M4 Max", "M2"]
multi_threaded = [29.1, 24.7, 72.1]
single_threaded = [383.1, 152.1, 186.9]

# Create positions for grouped bars
x_pos = 1:length(cpus)
bar_width = 0.35

# Create the plot with manual bar positioning
p = bar(x_pos .- bar_width/2, multi_threaded, 
        bar_width = bar_width,
        label = "Multi-threaded",
        color = :steelblue,
        title = "Solving the Arellano (2008) model",
        xlabel = "CPU",
        ylabel = "Time (seconds)",
        ylim   = (0,420),
        yticks = 0:60:420,
        xticks = (x_pos, cpus),
        legend = :topright,
        fontfamily = "Computer Modern",
        dpi = 300)

# Add single-threaded bars
bar!(x_pos .+ bar_width/2, single_threaded,
     bar_width = bar_width,
     label = "Single-threaded",
     color = :coral)

# Add value labels centered on the bars - using a simpler approach
for i in 1:length(cpus)
    # Multi-threaded bar labels
    annotate!(i - bar_width/2, multi_threaded[i]/2, 
             Plots.text("$(multi_threaded[i])s", 9,  :center, "Computer Modern"))
    # Single-threaded bar labels  
    annotate!(i + bar_width/2, single_threaded[i]/2, 
             Plots.text("$(single_threaded[i])s", 9, :center, "Computer Modern"))
end

display(p)

savefig(p, pwd()*"/files/benchmarks/arellano_2008.png")

#################################################
####### Processing large files #######
#################################################

# CPU performance data
cpus = ["Ryzen 9 9950X", "M4 Max", "M2"]
multi_threaded = [59.9, 79.0, 245.1]
single_threaded = [152.5, 346.4, 459.3]

# Create positions for grouped bars
x_pos = 1:length(cpus)
bar_width = 0.35

# Create the plot with manual bar positioning
p = bar(x_pos .- bar_width/2, multi_threaded, 
        bar_width = bar_width,
        label = "Multi-threaded",
        color = :steelblue,
        title = "Cleaning & processing a large dataset",
        xlabel = "CPU",
        ylabel = "Time (minutes)",
        ylim   = (0,480),
        yticks = 0:60:480,
        xticks = (x_pos, cpus),
        legend = :topleft,
        fontfamily = "Computer Modern",
        dpi = 300)

# Add single-threaded bars
bar!(x_pos .+ bar_width/2, single_threaded,
     bar_width = bar_width,
     label = "Single-threaded",
     color = :coral)

# Add value labels centered on the bars - using a simpler approach
for i in 1:length(cpus)
    # Multi-threaded bar labels
    annotate!(i - bar_width/2, multi_threaded[i]/2, 
             Plots.text("$(multi_threaded[i])m", 9,  :center, "Computer Modern"))
    # Single-threaded bar labels  
    annotate!(i + bar_width/2, single_threaded[i]/2, 
             Plots.text("$(single_threaded[i])m", 9, :center, "Computer Modern"))
end

display(p)

savefig(p, pwd()*"/files/benchmarks/lightcast.png")
