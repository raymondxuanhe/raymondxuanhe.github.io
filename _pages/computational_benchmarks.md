---
permalink: /computational_benchmarks/
title: "Computational Benchmarks"
author_profile: true
redirect_from: 
  - /computational_benchmarks/
  - /computational_benchmarks.html
---

Here, you can find the time to compute of three different CPUs for two computationally intensive tasks:

1. Reading, and cleaning 615 GB (compressed) of data,[^1] and 
2. Value function iteration solving the [Arellano (2008)](https://www.aeaweb.org/articles?id=10.1257/aer.98.3.690) model on a fine grid.[^2]

[comment]: #When I was deciding which computer I should buy (and if a PC, which parts), the vast majority of online reviews and benchmarks that I found made their comparisons and conclusions based on performance in video games, video encoding/creation and synthetic benchmarks (e.g. Cinebench and Geekbench). While these comparisons are fairly useful in terms of assessing components' *relative* performance, these tasks are very different to those that research economists typically undertake (especially if the task is multi-threaded), nor are they particularly informative of the real-world performance we can expect in our daily workflow. 

For full transparency, I have included all the relevant specificifications of each platform at the bottom of the page. All code was written on Julia. I am not sponsored by or affiliated with any company in any way.

## 1. Data analysis
![image info](/files/benchmarks/lightcast.png)

## 2. Value function iteration
![image info](/files/benchmarks/arellano_2008.png)

## Computer Specifications

### AMD Ryzen 9 9950X (16 cores, 32 threads - 30 threads used for parallelization)

* CPU Cooler: Arctic Liquid Freezer III 360 mm
* Motherboard: Gigabyte X870E Aorus Master
* GPU: Nvidia GeForce RTX 4080 Super (Founder's Edition)
* RAM: 64 GB of G.Skill DDR5 (2 x 32GB) 6000MT/s CL30
* SSD: 2 TB Samsung 990 Pro Internal SSD PCIe Gen 4x4 NVMe

### M4 Max (16 cores - 14 threads used for parallelization)

* 40-core GPU, 16-core neural engine
* 64 GB unified memory

### M2 Macbook Air (8 cores - 6 threads used for parallelization)

* 8-core GPU, 16-core neural engine
* 16 GB unified memory

[^1]: Parallelization was performed across individual files.
[^2]: The mesh consisted of 101 points for endowment and 951 points for bond holdings. Parallelization was performed across the endowment grid while evaluating the Bellman operator.
