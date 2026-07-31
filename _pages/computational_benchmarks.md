---
permalink: /computational_benchmarks/
title: "Computational Benchmarks"
author_profile: true
redirect_from: 
  - /computational_benchmarks.html
---

Here, I present the time-to-compute for solving the 
[Arellano (2008)](https://www.aeaweb.org/articles?id=10.1257/aer.98.3.690)
model via value function iteration on a fine grid on various processors.[^1] 
All code was written on Julia. I am not sponsored by or affiliated with any 
company in any way.

## CPUs
![image info](/files/benchmarks/cpu_benchmark.png)

## GPUs
![image info](/files/benchmarks/gpu_benchmark.png)

[^1]: The mesh consisted of 101 points for endowment and 951 points for debt. 
Parallelization was performed across the endowment grid and taking advantage of
Julia's column-major order.
