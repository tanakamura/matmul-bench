まとまるくん


# 結果
(理論値 : コア x SIMD x MulAdd x freq)

 * i7-4930K (3.40GHz, Ivy, 6core, 12thread, Linux)
  - SSE 版(!?) 1536x1536 106.37345[GFLOPS]  -mavx 付けないでビルドするほうがよい
  - 理論値 6 x 8 x 2 x 3.4 = 326.4[GFLOPS]
  - 32.6%
 * Zynq 7020 (667MHz, Cortex-A9, 2core, Linux)
  - A9 版 512x512 0.91347[GFLOPS]
  - 理論値 2 x 2 x 2 x 0.667 = 5.336[GFLOPS]
  - 17.1%
 * Celeron N2807 (1.58GHz, Silvermont, 2core, Linux)
  - SSE 版 512x512 7.76393[GFLOPS]
  - 理論値 2 x 4 x 2 x 1.58 = 25.28[GFLOPS]
  - 30.7%
 * Atom Z3740 (1.33GHz, Silvermont, 4core, Win 32bit)
  - SSE 版 512x512 11.09515[GFLOPS]
  - 理論値 4 x 4 x 2 x 1.33 = 42.56
  - 26.1%
 * Tegra4 (1.8GHz, Cortex-A15, 4core, Android)
  - NEON 版 512x512 8.85087[GFLOPS]
  - 理論値 4 x 4 x 2 x 1.8 = 57.6
  - 15.4%
 * i7-4770R (3.2GHz, Haswell, 4core, 8thread, Win 64bit)
  - FMA 版 1536x1536 104.93163[GFLOPS]
  - 理論値 4 x 8 x 4 x 3.2 = 409.6
  - 25.6%
 * Atom N550 (1.5GHz, Bonnell, 2core 4thred, Linux)
  - SSE 版 512x512 3.57132[GFLOPS]
  - 理論値 2 x 4 x 2 x 1.5 = 24.0
  - 14.9%
