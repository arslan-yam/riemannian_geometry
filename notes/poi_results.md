# POI setting: our POI index against the exact POI table

Queries are 1,000 random pairs among P random POIs, scored against exact MMP; single thread on the M1.
The table is built by one full VTP sweep per POI and stores P x P float32 distances, so it is
exact up to float32 rounding. Our index (eps 0.2, g 1 edge) builds the tree over the POIs
only; every answer was sound and none needed refining.

| mesh | n | P | table: build / size / query | ours: build / size (dense, 8-bit) / query / mean err | faster to build | smaller |
|---|---:|---:|---|---|---|---|
| HorseMount | 793 | 50 | 0.3 s / 0.01 MB / 0.033 µs | 1.5 s / 0.05, 0.04 MB / 0.32 µs / 0.068% | table | table |
| HorseMount | 793 | 200 | 1.1 s / 0.16 MB / 0.062 µs | 2.5 s / 0.17, 0.12 MB / 0.38 µs / 0.086% | table | ours |
| HorseMount | 793 | 793 | 4.0 s / 2.52 MB / 0.079 µs | 3.0 s / 1.00, 0.77 MB / 0.37 µs / 0.066% | ours | ours |
| RobinsonMountain | 3696 | 50 | 1.5 s / 0.01 MB / 0.027 µs | 14.7 s / 0.10, 0.08 MB / 0.27 µs / 0.057% | table | table |
| RobinsonMountain | 3696 | 200 | 5.3 s / 0.16 MB / 0.059 µs | 29.5 s / 0.32, 0.21 MB / 0.37 µs / 0.063% | table | table |
| RobinsonMountain | 3696 | 1000 | 31.4 s / 4.00 MB / 0.056 µs | 25.5 s / 1.49, 0.86 MB / 0.30 µs / 0.070% | ours | ours |
| RobinsonMountain | 3696 | 3000 | 78.6 s / 36.00 MB / 0.078 µs | 29.1 s / 5.53, 3.49 MB / 0.27 µs / 0.072% | ours | ours |
| Gunnison10k | 9610 | 50 | 2.9 s / 0.01 MB / 0.024 µs | 82.4 s / 0.15, 0.12 MB / 0.43 µs / 0.031% | table | table |
| Gunnison10k | 9610 | 200 | 17.9 s / 0.16 MB / 0.058 µs | 116.5 s / 0.47, 0.32 MB / 0.43 µs / 0.027% | table | table |
| Gunnison10k | 9610 | 1000 | 77.7 s / 4.00 MB / 0.041 µs | 117.2 s / 2.01, 1.07 MB / 0.35 µs / 0.021% | table | ours |
| Gunnison10k | 9610 | 3000 | 205.7 s / 36.00 MB / 0.092 µs | 149.4 s / 6.16, 3.11 MB / 0.30 µs / 0.028% | ours | ours |

Reading: below roughly 1,000 POIs on 3.7k–10k-vertex meshes the exact table wins on build time,
exactness and query time. Our index is cheaper to store from a few hundred POIs on, but its
build is dominated by the top cuts in both modes, so it is only cheaper to build once P is
large. The table grows as P^2 (2,500 POIs: 25 MB; 10,000: 400 MB), our index about linearly in P.
