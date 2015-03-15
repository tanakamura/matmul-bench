まとまるくん


# 結果
(理論値 : コア x SIMD x MulAdd x freq)

 * i7-4700MQ (2.4GHz, Haswell, 4C8T, Linux, TurboBoost なし)
  - FMA 版, 1920x1920, sec=0.05511, 256.85593[GFLOPS]
  - 理論値 4 x 8 x 4 x 2.4 = 307.2[GFLOPS]
  - 83.6%
