# All methods, one table per mesh

Single-threaded wall time on an Apple M1 (8 GB) for every method. `build` is index construction
(on-the-fly methods: none), `query` is per pair, errors are relative to exact MMP (|error| for
mean and max). Our C++ index and DGG/NeuroGF use the same 1,500 random vertex pairs per mesh (seed 7);
the EAR authors' harness uses 100 vertex pairs (seed 0, points 5 cm inside a face, float32 code);
on-the-fly exact solvers are timed on 5 full sweeps. Tags in brackets say whose code ran.
For the harness rows the logs give the signed mean, shown in absolute value (K-Algo and SE-Oracle
also answer below the truth, so their mean |error| is larger than shown).

## Summary: ours (C++, eps=0.1, g=0.5e, hybrid) against the best of each class

| mesh | n | ours: build / index / query / mean / max | EAR-Oracle sp=5 | DGG eps=1e-4 | exact on the fly (VTP) |
|---|---:|---|---|---|---|
| HorseMount | 793 | 3.1 s / 1.3 MB / 0.18 µs / 0.0235% / 1.16% | 23.5 s / 189.0 MB / 1.3 ms / 3.1016% / 24.51% | 730 ms / 4.4 MB / 1.4 ms / 0.0023% / 0.02% | 3.3 ms / exact |
| BigMount | 1440 | 10.4 s / 3.3 MB / 0.23 µs / 0.0220% / 0.86% | 52.4 s / 319.0 MB / 1.8 ms / 3.1914% / 6.41% | 1.5 s / 8.9 MB / 2.6 ms / 0.0025% / 0.02% | 6.6 ms / exact |
| HeadlightMount | 2470 | 29.1 s / 6.5 MB / 0.25 µs / 0.0156% / 1.30% | 107.3 s / 592.0 MB / 2.5 ms / 3.2905% / 8.47% | 2.9 s / 16.8 MB / 5.5 ms / 0.0026% / 0.02% | 12.0 ms / exact |
| RobinsonMountain | 3696 | 44.1 s / 10.7 MB / 0.28 µs / 0.0199% / 1.79% | 188.4 s / 672.0 MB / 3.8 ms / 3.5002% / 8.60% | 4.5 s / 25.3 MB / 7.9 ms / 0.0023% / 0.02% | 18.1 ms / exact |
| Gunnison10k | 9610 | 201.8 s / 37.4 MB / 0.36 µs / 0.0076% / 0.37% | 11.2 min / 702.0 MB / 10.2 ms / 3.4341% / 9.37% | 13.2 s / 72.8 MB / 22.3 ms / 0.0028% / 0.01% | 55.0 ms / exact |

Ours is certified on every query (LB <= d <= UB, sound on all pairs); none of the others
returns a per-query bound. Columns: build / index / query / mean |error| / max |error|.

## HorseMount

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| **ours, C++**, eps=0.2, g=1e, hybrid sweeps | 1.8 s | 1.00 MB | 0.20 µs | 0.0676% | 1.495% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 4.15% |
| **ours, C++**, eps=0.1, g=0.5e, hybrid sweeps | 3.1 s | 1.34 MB | 0.18 µs | 0.0235% | 1.159% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 3.41% |
| **ours, C++**, eps=0.1, g=0.5e, mmp sweeps | 3.8 s | 1.34 MB | 0.18 µs | 0.0235% | 1.159% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 3.41% |
| **ours, C++**, eps=0.1, g=0.5e, vtp sweeps | 3.8 s | 1.34 MB | 0.18 µs | 0.0235% | 1.159% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 3.41% |
| **ours, C++**, eps=0.05, g=0.25e, hybrid sweeps | 6.6 s | 2.01 MB | 0.22 µs | 0.0073% | 0.386% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 0.65% |
| ours, Python, eps=0.05 g=0.25e | 9.7 s | 2.01 MB | 1.23 µs | 0.0054% | 0.061% | 100 | certified, max UB err 0.14% |
| ours, Python, eps=0.1 g=0.5e | 4.6 s | 1.34 MB | 1.20 µs | 0.0114% | 0.135% | 100 | certified, max UB err 0.62% |
| ours, Python, eps=0.2 g=1e | 2.8 s | 1.00 MB | 1.22 µs | 0.0363% | 0.543% | 100 | certified, max UB err 1.67% |
| ours, A2A certified, eps=0.1 g=0.5e | 4.6 s | 1.34 MB | 122.90 µs | 0.1148% | 0.632% | 100 | A2A: sound True, max UB err 6.15%, width median 2.78% |
| EAR-Oracle, sp=10, eps=0.2 [EAR authors' C++] | 156.4 s | 541.00 MB | 5.7 ms | 2.5345% | 22.735% | 100 | no per-query bound |
| EAR-Oracle, sp=20, eps=0.2 [EAR authors' C++] | 484.0 s | 683.00 MB | 5.7 ms | 2.3239% | 18.611% | 100 | no per-query bound |
| EAR-Oracle, sp=5, eps=0.05 [EAR authors' C++] | 63.3 s | 348.00 MB | 8.3 ms | 3.1016% | 24.510% | 100 | no per-query bound |
| EAR-Oracle, sp=5, eps=0.2 [EAR authors' C++] | 23.5 s | 189.00 MB | 1.3 ms | 3.1016% | 24.510% | 100 | no per-query bound |
| K-Algo, eps=0.05 [EAR authors' C++] | — | — | 8.1 s | 0.0178% | 2.684% | 100 | no per-query bound |
| K-Algo, eps=0.2 [EAR authors' C++] | — | — | 8.1 s | 0.0178% | 2.684% | 100 | no per-query bound |
| MMP exact (on the fly, their CGAL code) [EAR authors' C++] | — | — | 34.4 ms | 0.0002% | 0.014% | 100 | no per-query bound |
| SE-Oracle, sp=5, eps=0.2 [EAR authors' C++] | 1.74 h | 942.00 MB | 11.6 ms | 0.8141% | 20.542% | 100 | no per-query bound |
| fixed Steiner (on the fly), sp=20 [EAR authors' C++] | — | — | 208.3 ms | 2.3121% | 18.611% | 100 | no per-query bound |
| fixed Steiner (on the fly), sp=5 [EAR authors' C++] | — | — | 26.4 ms | 3.1053% | 24.510% | 100 | no per-query bound |
| unfixed Steiner (on the fly), eps=0.05 [EAR authors' C++] | — | — | 23.2 s | 2.2827% | 19.918% | 100 | no per-query bound |
| unfixed Steiner (on the fly), eps=0.2 [EAR authors' C++] | — | — | 15.1 s | 2.3703% | 20.200% | 100 | no per-query bound |
| DGG, eps=0.0001 [authors' C++, macOS port] | 730 ms | 4.42 MB | 1.4 ms | 0.0023% | 0.018% | 1500 | graph + Dijkstra per query; eps bounds the mean error |
| DGG, eps=0.001 [authors' C++, macOS port] | 1.0 s | 1.76 MB | 724.00 µs | 0.0279% | 0.156% | 1500 | graph + Dijkstra per query; eps bounds the mean error |
| NeuroGF (distance branch, 30x100 iters) [authors' network] | 23.0 min | 0.25 MB | 280.24 µs | 0.4584% | 4.696% | 1500 | learned, no bound; batched 17.0 µs/pair |
| MMP exact, one sweep per query [Kirsanov, C++] | — | — | 9.0 ms | 0.0000% | 0.000% | 5 | exact, full sweep |
| VTP exact, one sweep per query [authors' C++] | — | — | 3.3 ms | 0.0000% | 0.000% | 5 | exact, full sweep |

## BigMount

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| **ours, C++**, eps=0.2, g=1e, hybrid sweeps | 6.1 s | 2.31 MB | 0.22 µs | 0.0667% | 1.708% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 6.20% |
| **ours, C++**, eps=0.1, g=0.5e, hybrid sweeps | 10.4 s | 3.26 MB | 0.23 µs | 0.0220% | 0.864% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 2.11% |
| **ours, C++**, eps=0.1, g=0.5e, mmp sweeps | 13.6 s | 3.26 MB | 0.23 µs | 0.0220% | 0.864% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 2.11% |
| **ours, C++**, eps=0.1, g=0.5e, vtp sweeps | 15.2 s | 3.26 MB | 0.23 µs | 0.0220% | 0.864% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 2.11% |
| **ours, C++**, eps=0.05, g=0.25e, hybrid sweeps | 22.2 s | 5.15 MB | 0.33 µs | 0.0067% | 0.292% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 1.13% |
| ours, Python, eps=0.05 g=0.25e | 32.1 s | 5.15 MB | 1.36 µs | 0.0049% | 0.064% | 100 | certified, max UB err 0.07% |
| ours, Python, eps=0.1 g=0.5e | 15.3 s | 3.26 MB | 1.26 µs | 0.0167% | 0.284% | 100 | certified, max UB err 0.28% |
| ours, Python, eps=0.2 g=1e | 8.9 s | 2.31 MB | 1.24 µs | 0.0788% | 1.826% | 100 | certified, max UB err 1.83% |
| ours, A2A certified, eps=0.1 g=0.5e | 15.3 s | 3.26 MB | 110.98 µs | 0.1442% | 2.687% | 100 | A2A: sound True, max UB err 20.20%, width median 3.53% |
| EAR-Oracle, sp=5, eps=0.2 [EAR authors' C++] | 52.4 s | 319.00 MB | 1.8 ms | 3.1914% | 6.412% | 100 | no per-query bound |
| K-Algo, eps=0.2 [EAR authors' C++] | — | — | 5.2 s | 0.0155% | 2.944% | 100 | no per-query bound |
| MMP exact (on the fly, their CGAL code) [EAR authors' C++] | — | — | 29.8 ms | 0.0002% | 0.004% | 100 | no per-query bound |
| fixed Steiner (on the fly), sp=5 [EAR authors' C++] | — | — | 21.2 ms | 3.1715% | 6.412% | 100 | no per-query bound |
| DGG, eps=0.0001 [authors' C++, macOS port] | 1.5 s | 8.89 MB | 2.6 ms | 0.0025% | 0.018% | 1500 | graph + Dijkstra per query; eps bounds the mean error |
| DGG, eps=0.001 [authors' C++, macOS port] | 460 ms | 3.42 MB | 1.3 ms | 0.0285% | 0.131% | 1500 | graph + Dijkstra per query; eps bounds the mean error |
| MMP exact, one sweep per query [Kirsanov, C++] | — | — | 18.8 ms | 0.0000% | 0.000% | 5 | exact, full sweep |
| VTP exact, one sweep per query [authors' C++] | — | — | 6.6 ms | 0.0000% | 0.000% | 5 | exact, full sweep |

## HeadlightMount

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| **ours, C++**, eps=0.2, g=1e, hybrid sweeps | 15.1 s | 4.43 MB | 0.22 µs | 0.0471% | 1.299% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 5.16% |
| **ours, C++**, eps=0.1, g=0.5e, hybrid sweeps | 29.1 s | 6.46 MB | 0.25 µs | 0.0156% | 1.300% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 2.82% |
| **ours, C++**, eps=0.1, g=0.5e, mmp sweeps | 35.8 s | 6.46 MB | 0.26 µs | 0.0156% | 1.300% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 2.82% |
| **ours, C++**, eps=0.1, g=0.5e, vtp sweeps | 29.5 s | 6.46 MB | 0.25 µs | 0.0156% | 1.300% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 2.82% |
| **ours, C++**, eps=0.05, g=0.25e, hybrid sweeps | 61.9 s | 10.52 MB | 0.40 µs | 0.0046% | 0.247% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 0.30% |
| ours, Python, eps=0.05 g=0.25e | 79.4 s | 10.52 MB | 1.48 µs | 0.0034% | 0.084% | 100 | certified, max UB err 0.09% |
| ours, Python, eps=0.1 g=0.5e | 37.1 s | 6.46 MB | 1.27 µs | 0.0170% | 0.357% | 100 | certified, max UB err 0.49% |
| ours, Python, eps=0.2 g=1e | 21.5 s | 4.43 MB | 1.29 µs | 0.0561% | 0.794% | 100 | certified, max UB err 1.23% |
| ours, A2A certified, eps=0.1 g=0.5e | 37.1 s | 6.46 MB | 88.91 µs | 0.0755% | 0.773% | 100 | A2A: sound True, max UB err 5.85%, width median 2.36% |
| EAR-Oracle, sp=5, eps=0.2 [EAR authors' C++] | 107.3 s | 592.00 MB | 2.5 ms | 3.2905% | 8.466% | 100 | no per-query bound |
| K-Algo, eps=0.2 [EAR authors' C++] | — | — | 8.9 s | 0.0123% | 0.036% | 20 | no per-query bound |
| MMP exact (on the fly, their CGAL code) [EAR authors' C++] | — | — | 72.9 ms | 0.0003% | 0.006% | 100 | no per-query bound |
| fixed Steiner (on the fly), sp=5 [EAR authors' C++] | — | — | 35.8 ms | 3.2749% | 8.466% | 100 | no per-query bound |
| DGG, eps=0.0001 [authors' C++, macOS port] | 2.9 s | 16.77 MB | 5.5 ms | 0.0026% | 0.016% | 1500 | graph + Dijkstra per query; eps bounds the mean error |
| DGG, eps=0.001 [authors' C++, macOS port] | 810 ms | 6.12 MB | 2.5 ms | 0.0285% | 0.132% | 1500 | graph + Dijkstra per query; eps bounds the mean error |
| MMP exact, one sweep per query [Kirsanov, C++] | — | — | 39.2 ms | 0.0000% | 0.000% | 5 | exact, full sweep |
| VTP exact, one sweep per query [authors' C++] | — | — | 12.0 ms | 0.0000% | 0.000% | 5 | exact, full sweep |

## RobinsonMountain

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| **ours, C++**, eps=0.2, g=1e, hybrid sweeps | 25.2 s | 7.33 MB | 0.23 µs | 0.0628% | 1.786% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 2.78% |
| **ours, C++**, eps=0.1, g=0.5e, hybrid sweeps | 44.1 s | 10.70 MB | 0.28 µs | 0.0199% | 1.786% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 1.79% |
| **ours, C++**, eps=0.1, g=0.5e, mmp sweeps | 54.7 s | 10.70 MB | 0.27 µs | 0.0199% | 1.786% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 1.79% |
| **ours, C++**, eps=0.1, g=0.5e, vtp sweeps | 50.8 s | 10.70 MB | 0.27 µs | 0.0199% | 1.786% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 1.79% |
| **ours, C++**, eps=0.05, g=0.25e, hybrid sweeps | 91.3 s | 17.47 MB | 0.43 µs | 0.0058% | 0.400% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 0.40% |
| ours, Python, eps=0.05 g=0.25e | 121.2 s | 17.47 MB | 1.48 µs | 0.0044% | 0.074% | 100 | certified, max UB err 0.07% |
| ours, Python, eps=0.1 g=0.5e | 59.0 s | 10.70 MB | 1.31 µs | 0.0232% | 0.384% | 100 | certified, max UB err 0.38% |
| ours, Python, eps=0.2 g=1e | 34.2 s | 7.33 MB | 1.30 µs | 0.0647% | 1.185% | 100 | certified, max UB err 2.21% |
| ours, A2A certified, eps=0.1 g=0.5e | 59.0 s | 10.70 MB | 84.95 µs | 0.1124% | 0.974% | 100 | A2A: sound True, max UB err 11.42%, width median 2.82% |
| EAR-Oracle, sp=5, eps=0.2 [EAR authors' C++] | 188.4 s | 672.00 MB | 3.8 ms | 3.5002% | 8.604% | 100 | no per-query bound |
| K-Algo, eps=0.2 [EAR authors' C++] | — | — | 20.8 s | 0.0967% | 2.229% | 20 | no per-query bound |
| MMP exact (on the fly, their CGAL code) [EAR authors' C++] | — | — | 108.5 ms | 0.0000% | 0.006% | 100 | no per-query bound |
| fixed Steiner (on the fly), sp=5 [EAR authors' C++] | — | — | 61.6 ms | 3.4842% | 8.604% | 100 | no per-query bound |
| DGG, eps=0.0001 [authors' C++, macOS port] | 4.5 s | 25.28 MB | 7.9 ms | 0.0023% | 0.018% | 1500 | graph + Dijkstra per query; eps bounds the mean error |
| DGG, eps=0.001 [authors' C++, macOS port] | 1.3 s | 9.35 MB | 3.8 ms | 0.0279% | 0.141% | 1500 | graph + Dijkstra per query; eps bounds the mean error |
| NeuroGF (distance branch, 30x100 iters) [authors' network] | 23.4 min | 0.25 MB | 260.12 µs | 0.6673% | 10.257% | 1500 | learned, no bound; batched 15.4 µs/pair |
| MMP exact, one sweep per query [Kirsanov, C++] | — | — | 57.5 ms | 0.0000% | 0.000% | 5 | exact, full sweep |
| VTP exact, one sweep per query [authors' C++] | — | — | 18.1 ms | 0.0000% | 0.000% | 5 | exact, full sweep |

## Gunnison10k

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| **ours, C++**, eps=0.2, g=1e, hybrid sweeps | 132.7 s | 23.97 MB | 0.24 µs | 0.0225% | 0.681% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 2.02% |
| **ours, C++**, eps=0.1, g=0.5e, hybrid sweeps | 201.8 s | 37.39 MB | 0.36 µs | 0.0076% | 0.374% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 0.94% |
| **ours, C++**, eps=0.1, g=0.5e, mmp sweeps | 283.7 s | 37.39 MB | 0.36 µs | 0.0076% | 0.374% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 0.94% |
| **ours, C++**, eps=0.1, g=0.5e, vtp sweeps | 215.0 s | 37.39 MB | 0.36 µs | 0.0076% | 0.374% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 0.94% |
| **ours, C++**, eps=0.05, g=0.25e, hybrid sweeps | 418.9 s | 64.22 MB | 0.61 µs | 0.0024% | 0.151% | 1500 | certified: LB≤d≤UB on all 1500 pairs (sound), max UB err 0.22% |
| ours, Python, eps=0.05 g=0.25e | 10.9 min | 64.22 MB | 1.64 µs | 0.0033% | 0.083% | 100 | certified, max UB err 0.08% |
| ours, Python, eps=0.1 g=0.5e | 320.1 s | 37.39 MB | 1.39 µs | 0.0112% | 0.301% | 100 | certified, max UB err 0.30% |
| ours, Python, eps=0.2 g=1e | 175.3 s | 23.97 MB | 1.28 µs | 0.0205% | 0.143% | 100 | certified, max UB err 0.20% |
| ours, A2A certified, eps=0.1 g=0.5e | 320.1 s | 37.39 MB | 78.28 µs | 0.0340% | 0.225% | 100 | A2A: sound True, max UB err 3.26%, width median 1.62% |
| EAR-Oracle, sp=5, eps=0.2 [EAR authors' C++] | 11.2 min | 702.00 MB | 10.2 ms | 3.4341% | 9.374% | 100 | no per-query bound |
| K-Algo, eps=0.2 [EAR authors' C++] | — | — | 80.7 s | 0.0134% | 0.074% | 20 | no per-query bound |
| MMP exact (on the fly, their CGAL code) [EAR authors' C++] | — | — | 478.7 ms | 0.0000% | 0.001% | 100 | no per-query bound |
| fixed Steiner (on the fly), sp=5 [EAR authors' C++] | — | — | 161.5 ms | 3.3714% | 9.374% | 100 | no per-query bound |
| DGG, eps=0.0001 [authors' C++, macOS port] | 13.2 s | 72.78 MB | 22.3 ms | 0.0028% | 0.013% | 1500 | graph + Dijkstra per query; eps bounds the mean error |
| DGG, eps=0.001 [authors' C++, macOS port] | 3.4 s | 25.30 MB | 10.7 ms | 0.0285% | 0.136% | 1500 | graph + Dijkstra per query; eps bounds the mean error |
| MMP exact, one sweep per query [Kirsanov, C++] | — | — | 225.7 ms | 0.0000% | 0.000% | 5 | exact, full sweep |
| VTP exact, one sweep per query [authors' C++] | — | — | 55.0 ms | 0.0000% | 0.000% | 5 | exact, full sweep |

## GunnisonForest100k

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| **ours, C++**, eps=0.2, g=1e, hybrid sweeps | 1.10 h | 533.76 MB | 0.51 µs | 0.0069% | 0.120% | 200 | certified: LB≤d≤UB on all 200 pairs (sound), max UB err 0.12% |
| ours, Python, eps=0.2 g=1e | 1.90 h | 533.75 MB | 1.56 µs | 0.0045% | 0.042% | 100 | certified, max UB err 0.07% |
| MMP exact (on the fly, their CGAL code) [EAR authors' C++] | — | — | 20.9 s | 0.0000% | 0.000% | 100 | no per-query bound |
| fixed Steiner (on the fly), sp=5 [EAR authors' C++] | — | — | 46.8 s | 2.9284% | 6.559% | 100 | no per-query bound |

## The UP-Oracle authors' oracles (POI setting)

Their oracles index a set of points of interest (POIs); queries are between POIs (200 random pairs).
Ours indexes every vertex and answers the same pairs. Their published A2A numbers are not
measured but extrapolated (POI results multiplied by vertices/POIs, `algorithms.h`); the
`every vertex a POI` rows are real runs.

### TJ (576 vertices), their 50 POIs

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| **ours, C++**, eps=0.2, g=1e (all vertices indexed) | 989 ms | 0.63 MB | 0.21 µs | 0.0890% | 2.181% | 200 | certified, sound=True, max UB err 2.18% |
| **ours, C++**, eps=0.1, g=0.5e (all vertices indexed) | 1.8 s | 0.88 MB | 0.16 µs | 0.0283% | 0.438% | 200 | certified, sound=True, max UB err 0.44% |
| **ours, C++**, eps=0.05, g=0.25e (all vertices indexed) | 4.7 s | 1.38 MB | 0.34 µs | 0.0067% | 0.107% | 200 | certified, sound=True, max UB err 0.16% |
| EAR-Oracle [UP-Oracle authors' C++] | 2.0 s | 26.45 MB | 8.7 ms | 0.0118% | 0.109% | 200 | POI oracle |
| ESP-Fly [UP-Oracle authors' C++] | 0 ms | — | 8.1 ms | 0.3476% | 0.879% | 200 | on the fly |
| RC-Oracle [UP-Oracle authors' C++] | 247 ms | 0.97 MB | 1.70 µs | 0.0000% | 0.000% | 200 | POI oracle |
| UP-Oracle [UP-Oracle authors' C++] | 276 ms | 0.06 MB | 3.10 µs | 2.3049% | 13.415% | 200 | POI oracle |
| WAV-Fly [UP-Oracle authors' C++] | 0 ms | — | 5.1 ms | 0.0000% | 0.000% | 200 | on the fly |
| WSPD-Oracle [UP-Oracle authors' C++] | 271 ms | 1.08 MB | 1.70 µs | 0.0000% | 0.000% | 200 | POI oracle |

### TJ, every vertex a POI (576)

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| **ours, C++**, eps=0.2, g=1e (all vertices indexed) | 1.4 s | 0.63 MB | 0.19 µs | 0.0886% | 1.417% | 200 | certified, sound=True, max UB err 1.45% |
| **ours, C++**, eps=0.1, g=0.5e (all vertices indexed) | 1.8 s | 0.88 MB | 0.18 µs | 0.0284% | 0.465% | 200 | certified, sound=True, max UB err 0.46% |
| **ours, C++**, eps=0.05, g=0.25e (all vertices indexed) | 4.0 s | 1.38 MB | 0.23 µs | 0.0070% | 0.162% | 200 | certified, sound=True, max UB err 0.16% |
| EAR-Oracle [UP-Oracle authors' C++] | 5.6 s | 48.06 MB | 12.7 ms | 0.0371% | 2.285% | 200 | POI oracle |
| ESP-Fly [UP-Oracle authors' C++] | 0 ms | — | 8.3 ms | 0.2962% | 1.023% | 200 | on the fly |
| RC-Oracle [UP-Oracle authors' C++] | 4.7 s | 112.97 MB | 6.30 µs | 0.1427% | 8.695% | 200 | POI oracle |
| UP-Oracle [UP-Oracle authors' C++] | 15.0 s | 9.85 MB | 292.90 µs | 1.3410% | 7.856% | 200 | POI oracle |
| WAV-Fly [UP-Oracle authors' C++] | 0 ms | — | 4.8 ms | 0.0000% | 0.000% | 200 | on the fly |
| WSPD-Oracle [UP-Oracle authors' C++] | 5.1 s | 125.52 MB | 5.20 µs | 0.1427% | 8.695% | 200 | POI oracle |

### HorseMount, 50 random POIs

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| **ours, C++**, eps=0.2, g=1e (all vertices indexed) | 1.8 s | 1.00 MB | 0.16 µs | 0.0771% | 0.802% | 200 | certified, sound=True, max UB err 2.74% |
| **ours, C++**, eps=0.1, g=0.5e (all vertices indexed) | 3.1 s | 1.34 MB | 0.16 µs | 0.0252% | 0.447% | 200 | certified, sound=True, max UB err 1.29% |
| **ours, C++**, eps=0.05, g=0.25e (all vertices indexed) | 6.6 s | 2.01 MB | 0.21 µs | 0.0070% | 0.121% | 200 | certified, sound=True, max UB err 0.13% |
| EAR-Oracle [UP-Oracle authors' C++] | 1.8 s | 8.36 MB | 8.1 ms | 0.0394% | 1.275% | 200 | POI oracle |
| ESP-Fly [UP-Oracle authors' C++] | 0 ms | — | 12.7 ms | 0.4192% | 2.411% | 200 | on the fly |
| RC-Oracle [UP-Oracle authors' C++] | 430 ms | 1.00 MB | 2.40 µs | 0.0317% | 3.167% | 200 | POI oracle |
| UP-Oracle [UP-Oracle authors' C++] | 493 ms | 0.12 MB | 4.20 µs | 2.5266% | 9.949% | 200 | POI oracle |
| WAV-Fly [UP-Oracle authors' C++] | 0 ms | — | 9.1 ms | 0.0000% | 0.000% | 200 | on the fly |
| WSPD-Oracle [UP-Oracle authors' C++] | 473 ms | 1.11 MB | 3.00 µs | 0.0317% | 3.167% | 200 | POI oracle |

### HorseMount, every vertex a POI (793)

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| **ours, C++**, eps=0.2, g=1e (all vertices indexed) | 1.8 s | 1.00 MB | 0.19 µs | 0.0564% | 0.562% | 200 | certified, sound=True, max UB err 5.95% |
| **ours, C++**, eps=0.1, g=0.5e (all vertices indexed) | 3.2 s | 1.34 MB | 0.18 µs | 0.0226% | 0.535% | 200 | certified, sound=True, max UB err 1.03% |
| **ours, C++**, eps=0.05, g=0.25e (all vertices indexed) | 6.6 s | 2.01 MB | 0.22 µs | 0.0103% | 0.535% | 200 | certified, sound=True, max UB err 0.73% |
| EAR-Oracle [UP-Oracle authors' C++] | 11.4 s | 61.37 MB | 8.2 ms | 0.0342% | 0.401% | 200 | POI oracle |
| ESP-Fly [UP-Oracle authors' C++] | 0 ms | — | 12.5 ms | 0.4087% | 1.137% | 200 | on the fly |
| RC-Oracle [UP-Oracle authors' C++] | 11.8 s | 219.74 MB | 7.20 µs | 0.5948% | 10.208% | 200 | POI oracle |
| UP-Oracle [UP-Oracle authors' C++] | 26.7 s | 27.77 MB | 748.80 µs | 0.8027% | 5.843% | 200 | POI oracle |
| WAV-Fly [UP-Oracle authors' C++] | 0 ms | — | 9.2 ms | 0.0000% | 0.000% | 200 | on the fly |
| WSPD-Oracle [UP-Oracle authors' C++] | 13.0 s | 244.16 MB | 6.80 µs | 0.5948% | 10.208% | 200 | POI oracle |

## Weighted HorseMount (face costs in (1, 5)), reference = best valid path any method found

| method | build | index | query | mean err | max err | pairs | notes |
|---|---:|---:|---:|---:|---:|---:|---|
| ours: Steiner graph m=24 Dijkstra (no index) | — | — | 99.3 ms | 0.0027% | 0.257% | 100 | weighted; UB is a path cost |
| ours: certified eps=0.1 g=0.5e, Steiner m=10 (UB) | 32.8 s | 2.90 MB | 1.31 µs | 0.2299% | 2.083% | 100 | weighted; UB is a path cost |
| ours: certified eps=0.1 g=0.5e, Steiner m=10 (estimate) | 32.8 s | 2.90 MB | 1.32 µs | 0.1960% | 2.083% | 100 | weighted; UB is a path cost |
| ours: certified eps=0.1 g=0.5e, Steiner m=5 (UB) | 10.9 s | 2.90 MB | 1.32 µs | 0.3856% | 3.107% | 100 | weighted; UB is a path cost |
| ours: certified eps=0.1 g=0.5e, Steiner m=5 (estimate) | 10.9 s | 2.90 MB | 1.31 µs | 0.3459% | 3.107% | 100 | weighted; UB is a path cost |
| ours: certified eps=0.2 g=1e, Steiner m=10 (UB) | 18.7 s | 1.79 MB | 1.22 µs | 0.5943% | 4.489% | 100 | weighted; UB is a path cost |
| ours: certified eps=0.2 g=1e, Steiner m=10 (estimate) | 18.7 s | 1.79 MB | 1.22 µs | 0.5458% | 4.489% | 100 | weighted; UB is a path cost |
| ours: certified eps=0.2 g=1e, Steiner m=5 (UB) | 6.6 s | 1.79 MB | 1.23 µs | 0.7935% | 4.703% | 100 | weighted; UB is a path cost |
| ours: certified eps=0.2 g=1e, Steiner m=5 (estimate) | 6.6 s | 1.79 MB | 1.24 µs | 0.7395% | 4.703% | 100 | weighted; UB is a path cost |
| ours: dijkstra on mesh graph (no index) | — | — | 783.21 µs | 3.3631% | 48.532% | 100 | weighted; UB is a path cost |
| ours: h2h | 78 ms | 0.65 MB | 5.65 µs | 3.3631% | 48.532% | 100 | weighted; UB is a path cost |
| EAR-Oracle, sp=10, eps=0.2 [EAR authors' C++] | 91.1 s | 299.00 MB | 1.6 ms | 2.2721% | 48.584% | 100 | no per-query bound |
| EAR-Oracle, sp=20, eps=0.2 [EAR authors' C++] | 305.3 s | 398.00 MB | 3.5 ms | 2.1709% | 48.584% | 100 | no per-query bound |
| EAR-Oracle, sp=5, eps=0.2 [EAR authors' C++] | 25.9 s | 166.00 MB | 1.1 ms | 2.4239% | 48.584% | 100 | no per-query bound |
| fixed Steiner (on the fly), sp=20 [EAR authors' C++] | — | — | 142.4 ms | 2.1600% | 48.584% | 100 | no per-query bound |
| fixed Steiner (on the fly), sp=5 [EAR authors' C++] | — | — | 11.8 ms | 2.4264% | 48.584% | 100 | no per-query bound |
| unfixed Steiner (on the fly), eps=0.05 [EAR authors' C++] | — | — | 10.9 s | 1.8567% | 5.225% | 30 | no per-query bound |
