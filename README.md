まとまるくん


# 結果
(理論値 : コア x SIMD x MulAdd x freq)

 * i7-4930K (3.40GHz, Ivy, 6core, 12thread, Linux)
  - SSE 版(!?) 1536x1536 106.37345[GFLOPS]  -mavx 付けないでビルドするほうがよい
  - 理論値 6 x 8 x 2 x 3.4 = 326.4G[FLOPS]
  - 32.6%
 * Cortex A9 (667MHz, A9, 2core, Linux)
  - A9 版 512x512 0.91347[GFLOPS]
  - 理論値 2 x 2 x 2 x 0.667 = 5.336[GFLOPS]
  - 17.1%
 * Celeron N2807 (1.58GHz, Silvermont, 2core, Linux)
  - SSE 版 512x512 7.76393[GFLOPS]
  - 理論値 2 x 4 x 2 x 1.58 = 25.28[GFLOPS]
  - 30.7%
 * Atom Z3740 (1.33GHz, Silvermont, 4core, Win)
  - SSE 版 512x512 11.09515[GFLOPS]
  - 理論値 4 x 4 x 2 x 1.33 = 42.56
  - 26.1%
