まとまるくん


# 結果
(理論値 : コア x SIMD x MulAdd x freq)

 * i7-4930K (3.40GHz, Ivy, 6core, 12thread, Linux, Turbo 3.9GHz あり)
  - AVX 版 1536x1536 174.42793[GFLOPS]
  - 理論値 6 x 8 x 2 x 3.4 = 326.4[GFLOPS]
  - -std=gnu99 -Wall -O2 -fopenmp -ffast-math -mtune=native -save-temps -march=native
  - 53.4%
 * Zynq 7020 (667MHz, Cortex-A9, 2core, Linux)
  - NEON 版 512x512 1.25446[GFLOPS]
  - 理論値 2 x 2 x 2 x 0.667 = 5.336[GFLOPS]
  - -std=gnu99 -Wall -O2 -ffast-math -mfloat-abi=hard -fopenmp -mfpu=neon -mtune=cortex-a9 -mcpu=cortex-a9 -falign-loops=16 
  - 23.5%
 * Celeron N2807 (1.58GHz, Silvermont, 2core, Linux)
  - SSE 版 512x512 7.76393[GFLOPS]
  - 理論値 2 x 4 x 2 x 1.58 = 25.28[GFLOPS]
  - 30.7%
 * Atom Z3740 (1.33GHz, Silvermont, 4core, Win 32bit)
  - SSE 版 512x512 11.09515[GFLOPS]
  - 理論値 4 x 4 x 2 x 1.33 = 42.56
  - 26.1%
 * Tegra4 (1.8GHz, Cortex-A15, 4core, Android) (なんか4コア使うとクロックは1.6GHzになる)
  - NEON 版 1024x1024 9.80157[GFLOPS]
  - 理論値 4 x 4 x 2 x 1.8 = 57.6
  - -std=gnu99 -Wall -O2 -ffast-math -mfloat-abi=softfp -mfpu=neon-vfpv4 -mtune=cortex-a15 -mcpu=cortex-a15 -fopenmp
  - 17.0%
 * i7-4770R (3.2GHz, Haswell, 4core, 8thread, Linux, Turbo 3.9GHz あり)
  - FMA 版 1536x1536 203.07178[GFLOPS]
  - 理論値 4 x 8 x 4 x 3.2 = 409.6
  - 49.6%
 * Atom N550 (1.5GHz, Bonnell, 2core 4thred, Linux)
  - SSE 版 512x512 3.57132[GFLOPS]
  - 理論値 2 x 4 x 2 x 1.5 = 24.0
  - 14.9%
