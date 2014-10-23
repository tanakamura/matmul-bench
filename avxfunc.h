static void
AVX_FUNC_NAME(float * __restrict out,
              const float* __restrict inL,
              const float* __restrict inR,
              unsigned int n)
{
    unsigned long block_size = 16;
    unsigned long block_size_k = 64;

    int i00;

#pragma omp parallel for schedule(dynamic)
    for (i00=0; i00<n; i00+=block_size) {
        for (int j0=0; j0<n; j0+=block_size) {
            for (int k0=0; k0<n; k0+=block_size_k) {
                for (int bi=0; bi<block_size; bi+=4) {
                    int i0 = i00+bi+0;
                    int i1 = i00+bi+1;
                    int i2 = i00+bi+2;
                    int i3 = i00+bi+3;

                    float *outp0 = &out[i0*n+j0];
                    float *outp1 = &out[i1*n+j0];
                    float *outp2 = &out[i2*n+j0];
                    float *outp3 = &out[i3*n+j0];

                    __m256 *outp8_0 = (__m256*)outp0;
                    __m256 *outp8_1 = (__m256*)outp1;
                    __m256 *outp8_2 = (__m256*)outp2;
                    __m256 *outp8_3 = (__m256*)outp3;

                    __m256 vout0_0 = _mm256_setzero_ps();
                    __m256 vout0_1 = _mm256_setzero_ps();
                    __m256 vout1_0 = _mm256_setzero_ps();
                    __m256 vout1_1 = _mm256_setzero_ps();
                    __m256 vout2_0 = _mm256_setzero_ps();
                    __m256 vout2_1 = _mm256_setzero_ps();
                    __m256 vout3_0 = _mm256_setzero_ps();
                    __m256 vout3_1 = _mm256_setzero_ps();

                    for (long bk=0; bk<block_size_k; bk+=8) {
                        //float lik0, lik1, lik2, lik3;

#define AVX_K8(K)                                                       \
                        {                                               \
                            long k = k0+bk+K;                           \
                            const float *inRp = &inR[k*n+j0];           \
                                                                        \
                            __m256 vr0 = _mm256_load_ps(&inRp[0]);      \
                            __m256 vr1 = _mm256_load_ps(&inRp[8]);      \
                                                                        \
                            _mm_prefetch((const char*)(inRp + n*3), _MM_HINT_T0); \
                                                                        \
                            AVX_OP(0,0);                                \
                            AVX_OP(0,1);                                \
                                                                        \
                            AVX_OP(1,0);                                \
                            AVX_OP(1,1);                                \
                                                                        \
                            AVX_OP(2,0);                                \
                            AVX_OP(2,1);                                \
                                                                        \
                            AVX_OP(3,0);                                \
                            AVX_OP(3,1);                                \
                        }

                        //__m256 vlik0 = *(__m256*)&inL[i0*n+k0+bk];
                        //__m256 vlik1 = *(__m256*)&inL[i1*n+k0+bk];
                        //__m256 vlik2 = *(__m256*)&inL[i2*n+k0+bk];
                        //__m256 vlik3 = *(__m256*)&inL[i3*n+k0+bk];

                        __m256 lik0_8, lik1_8, lik2_8, lik3_8;

#if 0
#define AVX_K8_L(K)                             \
                        lik0_8 = _mm256_shuffle_ps(vlik0,vlik0,_MM_SHUFFLE(K,K,K,K)); \
                        lik1_8 = _mm256_shuffle_ps(vlik1,vlik1,_MM_SHUFFLE(K,K,K,K)); \
                        lik2_8 = _mm256_shuffle_ps(vlik2,vlik2,_MM_SHUFFLE(K,K,K,K)); \
                        lik3_8 = _mm256_shuffle_ps(vlik3,vlik3,_MM_SHUFFLE(K,K,K,K)); \

#else

#define AVX_K8_L(K)                             \
                        lik0_8 = _mm256_set1_ps(inL[i0*n+k0+bk+K]); \
                        lik1_8 = _mm256_set1_ps(inL[i1*n+k0+bk+K]); \
                        lik2_8 = _mm256_set1_ps(inL[i2*n+k0+bk+K]); \
                        lik3_8 = _mm256_set1_ps(inL[i3*n+k0+bk+K]); \

#endif

                        AVX_K8_L(0);
                        AVX_K8(0);
                        AVX_K8_L(1);
                        AVX_K8(1);
                        AVX_K8_L(2);
                        AVX_K8(2);
                        AVX_K8_L(3);
                        AVX_K8(3);

                        //vlik0 = _mm256_permute2f128_ps(vlik0,vlik0,0x11);
                        //vlik1 = _mm256_permute2f128_ps(vlik1,vlik1,0x11);
                        //vlik2 = _mm256_permute2f128_ps(vlik2,vlik2,0x11);
                        //vlik3 = _mm256_permute2f128_ps(vlik3,vlik3,0x11);

                        AVX_K8_L(4);
                        AVX_K8(4);
                        AVX_K8_L(5);
                        AVX_K8(5);
                        AVX_K8_L(6);
                        AVX_K8(6);
                        AVX_K8_L(7);
                        AVX_K8(7);
                    }

                    if (k0 == 0) {
                        outp8_0[0] = vout0_0;
                        outp8_0[1] = vout0_1;
                        outp8_1[0] = vout1_0;
                        outp8_1[1] = vout1_1;
                        outp8_2[0] = vout2_0;
                        outp8_2[1] = vout2_1;
                        outp8_3[0] = vout3_0;
                        outp8_3[1] = vout3_1;
                    } else {
                        outp8_0[0] = _mm256_add_ps(outp8_0[0], vout0_0);
                        outp8_0[1] = _mm256_add_ps(outp8_0[1], vout0_1);
                        outp8_1[0] = _mm256_add_ps(outp8_1[0], vout1_0);
                        outp8_1[1] = _mm256_add_ps(outp8_1[1], vout1_1);
                        outp8_2[0] = _mm256_add_ps(outp8_2[0], vout2_0);
                        outp8_2[1] = _mm256_add_ps(outp8_2[1], vout2_1);
                        outp8_3[0] = _mm256_add_ps(outp8_3[0], vout3_0);
                        outp8_3[1] = _mm256_add_ps(outp8_3[1], vout3_1);
                    }
                }
            }
        }
    }
}

#undef AVX_OP
#undef AVX_FUNC_NAME

