static void
AVX_FUNC_NAME(float * __restrict out,
              const float* __restrict inL,
              const float* __restrict inR,
              unsigned int n)
{
    unsigned int block_size = 32;
    int i00, i;

#pragma omp parallel for schedule(dynamic)
    for (i00=0; i00<n; i00+=block_size) {
        for (int j0=0; j0<n; j0+=block_size) {
            for (int k0=0; k0<n; k0+=block_size) {
                for (int bi=0; bi<block_size; bi+=2) {
                    int i0 = i00+bi+0;
                    int i1 = i00+bi+1;

                    float *outp0 = &out[i0*n+j0];
                    float *outp1 = &out[i1*n+j0];

                    __m256 *outp8_0 = (__m256*)outp0;
                    __m256 *outp8_1 = (__m256*)outp1;

                    __m256 vout0_0 = _mm256_setzero_ps();
                    __m256 vout0_1 = _mm256_setzero_ps();
                    __m256 vout0_2 = _mm256_setzero_ps();
                    __m256 vout0_3 = _mm256_setzero_ps();

                    __m256 vout1_0 = _mm256_setzero_ps();
                    __m256 vout1_1 = _mm256_setzero_ps();
                    __m256 vout1_2 = _mm256_setzero_ps();
                    __m256 vout1_3 = _mm256_setzero_ps();

                    if (k0 == 0) {
                        _mm_prefetch((const char*)(outp8_0) ,_MM_HINT_T0);
                        _mm_prefetch((const char*)(outp8_1) ,_MM_HINT_T0);
                    }

                    for (long bk=0; bk<block_size; bk++) {
                        long k = k0+bk;

                        const float *inRp = &inR[k*n+j0];

                        _mm_prefetch((const char*)(inRp + n), _MM_HINT_T0);

                        float lik0 = inL[i0*n+k];
                        float lik1 = inL[i1*n+k];

                        __m256 lik0_8 = _mm256_set1_ps(lik0);
                        __m256 lik1_8 = _mm256_set1_ps(lik1);


                        AVX_OP(0,0);
                        AVX_OP(0,1);
                        AVX_OP(0,2);
                        AVX_OP(0,3);

                        AVX_OP(1,0);
                        AVX_OP(1,1);
                        AVX_OP(1,2);
                        AVX_OP(1,3);
                    }

                    if (k0 == 0) {
                        outp8_0[0] = vout0_0;
                        outp8_0[1] = vout0_1;
                        outp8_0[2] = vout0_2;
                        outp8_0[3] = vout0_3;

                        outp8_1[0] = vout1_0;
                        outp8_1[1] = vout1_1;
                        outp8_1[2] = vout1_2;
                        outp8_1[3] = vout1_3;
                    } else {
                        outp8_0[0] = _mm256_add_ps(outp8_0[0], vout0_0);
                        outp8_0[1] = _mm256_add_ps(outp8_0[1], vout0_1);
                        outp8_0[2] = _mm256_add_ps(outp8_0[2], vout0_2);
                        outp8_0[3] = _mm256_add_ps(outp8_0[3], vout0_3);

                        outp8_1[0] = _mm256_add_ps(outp8_1[0], vout1_0);
                        outp8_1[1] = _mm256_add_ps(outp8_1[1], vout1_1);
                        outp8_1[2] = _mm256_add_ps(outp8_1[2], vout1_2);
                        outp8_1[3] = _mm256_add_ps(outp8_1[3], vout1_3);
                    }
                }
            }
        }
    }
}

#undef AVX_OP
#undef AVX_FUNC_NAME

