---
permalink: /computational_benchmarks/
title: "Computational Benchmarks"
author_profile: true
redirect_from: 
  - /computational_benchmarks/
  - /computational_benchmarks.html
---

Here, you can find the time to compute of three different CPUs for two computationally intensive tasks:

1. Reading, and cleaning 900GB (compressed) of data[^1] 
2. Value function iteration solving the [Arellano (2008)](https://www.aeaweb.org/articles?id=10.1257/aer.98.3.690) model on a fine grid[^2]

[comment]: #When I was deciding which computer I should buy (and if a PC, which parts), the vast majority of online reviews and benchmarks that I found made their comparisons and conclusions based on performance in video games, video encoding/creation and synthetic benchmarks (e.g. Cinebench and Geekbench). While these comparisons are fairly useful in terms of assessing components' *relative* performance, these tasks are very different to those that research economists typically undertake (especially if the task is multi-threaded), nor are they particularly informative of the real-world performance we can expect in our daily workflow. 

For full transparency, I have included all the relevant specificifications of each platform at the bottom of the page. All code was written on Julia. I am not sponsored by or affiliated with any company in any way.

## 1. Data analysis




## 2. Value function iteration



## Computer Specifications

### AMD Ryzen 9 9950X (16 cores, 32 threads)

* CPU Cooler: Arctic Liquid Freezer III 360 mm
* Motherboard: Gigabyte X870E Aorus Master
* GPU: Nvidia GeForce RTX 4080 Super (Founder's Edition)
* RAM: 64 GB of G.Skill DDR5 (2 x 32GB) 6000MT/s CL30
* SSD: 2TB Samsung 990 Pro Internal SSD PCIe Gen 4x4 NVMe

### AMD Ryzen 9 5900X (12 cores, 24 threads)

* CPU Cooler: Corsair iCUE H100i ELITE CAPELLIX XT 240mm
* Motherboard: MSI MAG B550 Tomahawk
* GPU: Nvidia GeForce RTX 3060 (Gigabyte)
* RAM: 64GB of DDR4 Corsair Vengeance Pro 3600MT/s CL18
* SSD: 2TB Corsair MP600 Elite PCIe Gen 4x4 NVMe

### M2 Macbook Air (8 cores)

* 8-core GPU, 16-core neural engine
* 16GB unified memory

[^1]: Parallelization was performed across individual files.
[^2]: The mesh consisted of 101 points for endowment and 951 points for bond holdings. Parallelization was performed across the endowment grid while evaluating the Bellman operator.
