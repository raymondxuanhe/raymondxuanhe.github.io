---
permalink: /computational_benchmarks/
title: "Computational Benchmarks"
author_profile: true
redirect_from: 
  - /computational_benchmarks.html
---

I present the various processors' time-to-compute for solving the
[Arellano (2008)](https://www.aeaweb.org/articles?id=10.1257/aer.98.3.690)
model via value function iteration on a discrete grid.[^1] All code was 
written on Julia. I am not sponsored by or affiliated with any company 
in any way.

[comment]: #When I was deciding which computer I should buy (and if a PC, 
which parts), the vast majority of online reviews and benchmarks that I 
found made their comparisons and conclusions based on performance in video 
games, video encoding/creation and synthetic benchmarks (e.g. Cinebench and 
Geekbench). While these comparisons are fairly useful in terms of assessing 
components' *relative* performance, these tasks are very different to those t
hat research economists typically undertake (especially if the task is 
multi-threaded), nor are they particularly informative of the real-world 
performance we can expect in our daily workflow.

## CPUs
![image info](/files/benchmarks/cpu_benchmark.png)

## GPUs
![image info](/files/benchmarks/gpu_benchmark.png)

[^1]: The mesh consisted of 101 points for endowment and 951 points for debt. 
Parallelization was performed across the endowment grid and taking advantage of
Julia's column-major order.
