---
permalink: /computational_benchmarks/
title: "Computational Benchmarks"
author_profile: true
redirect_from: 
  - /computational_benchmarks.html
---

You can find the time-to-compute of four different processors for solving the
[Arellano (2008)](https://www.aeaweb.org/articles?id=10.1257/aer.98.3.690)
via value function iteration.[^1] All code was written on Julia. I am not sponsored by or affiliated with any company in any way.

[comment]: #When I was deciding which computer I should buy (and if a PC, which parts), the vast majority of online reviews and benchmarks that I found made their comparisons and conclusions based on performance in video games, video encoding/creation and synthetic benchmarks (e.g. Cinebench and Geekbench). While these comparisons are fairly useful in terms of assessing components' *relative* performance, these tasks are very different to those that research economists typically undertake (especially if the task is multi-threaded), nor are they particularly informative of the real-world performance we can expect in our daily workflow. 

[comment]: #For full transparency, I have included all the relevant specificifications of each platform at the bottom of the page. All code was written on Julia. I am not sponsored by or affiliated with any company in any way.

## CPUs
![image info](/files/benchmarks/cpu_benchmarks.png)

## GPUs
![image info](/files/benchmarks/cpu_benchmarks.png)

[^1]: The mesh consisted of 101 points for endowment and 951 points for debt. Parallelization was performed across the endowment grid.