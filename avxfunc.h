#ifndef _mm256_set_m128
#define _mm256_set_m128(hi,lo) \
    _mm256_insertf128_ps(_mm256_castps128_ps256(lo), (hi), 0x1)
#endif


static void
CONCAT(AVX_FUNC_NAME,_)(long i00,
                        long j0,
                        long k0,
                        long bi,
                        float * __restrict out,
                        const float* __restrict inL,
                        const float* __restrict inR,
                        long n)
{
    long i0 = i00+bi+0;
    long i1 = i00+bi+1;
    long i2 = i00+bi+2;
    long i3 = i00+bi+3;

    float *outp0 = &out[i0*n+j0];
    float *outp1 = &out[i1*n+j0];
    float *outp2 = &out[i2*n+j0];
    float *outp3 = &out[i3*n+j0];

    _mm_prefetch((const char*)(outp0), _MM_HINT_T0);
    _mm_prefetch((const char*)(outp1), _MM_HINT_T0);
    _mm_prefetch((const char*)(outp2), _MM_HINT_T0);
    _mm_prefetch((const char*)(outp3), _MM_HINT_T0);

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


#define AVX_K8(K)                                               \
    {                                                           \
        __m256 vr0 = _mm256_load_ps(&inR0[0]);                  \
        __m256 vr1 = _mm256_load_ps(&inR0[8]);                  \
                                                                \
        _mm_prefetch((const char*)(inR0 + n*4), _MM_HINT_T0);   \
        inR0 += n;                                              \
                                                                \
        AVX_OP(0,0);                                            \
        AVX_OP(0,1);                                            \
                                                                \
        AVX_OP(1,0);                                            \
        AVX_OP(1,1);                                            \
                                                                \
        AVX_OP(2,0);                                            \
        AVX_OP(2,1);                                            \
                                                                \
        AVX_OP(3,0);                                            \
        AVX_OP(3,1);                                            \
    }


#define AVX_K8_L(K)                             \
    lik0_8 = _mm256_set1_ps(inL[i0*n+k0+bk+K]); \
    lik1_8 = _mm256_set1_ps(inL[i1*n+k0+bk+K]); \
    lik2_8 = _mm256_set1_ps(inL[i2*n+k0+bk+K]); \
    lik3_8 = _mm256_set1_ps(inL[i3*n+k0+bk+K]);


#define AVX_K8_L2(K)                            \
    lik0_8 = _mm256_set1_ps(inL0[n*0]);         \
    lik1_8 = _mm256_set1_ps(inL0[n*1]);         \
    lik2_8 = _mm256_set1_ps(inL0[n*2]);         \
    lik3_8 = _mm256_set1_ps(inL0[n*3]);         \
    inL0++;

#if 1
    __m256 lik0_8, lik1_8, lik2_8, lik3_8;

#define LC(N) AVX_K8_L2(N); AVX_K8(N);

    const float *inL0 = &inL[i0*n+k0];
    const float *inR0 = &inR[k0*n+j0];

    M16(LC);
    M16(LC);
    M16(LC);
    M16(LC);

#else
    for (long bk=0; bk<block_size_k; bk++) {
        __m256 lik0_8, lik1_8, lik2_8, lik3_8;

        AVX_K8_L_LO(0);
        AVX_K8(0);
    }
#endif

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


static void
AVX_FUNC_NAME(float * __restrict out,
              const float* __restrict inL,
              const float* __restrict inR,
              unsigned int n)
{
    unsigned long block_size = 16;
    unsigned long block_size_k = 64;

    long i00;

#pragma omp parallel for schedule(dynamic)
    for (i00=0; i00<n; i00+=block_size) {
        for (long j0=0; j0<n; j0+=block_size) {
            for (long k0=0; k0<n; k0+=block_size_k) {
                for (long bi=0; bi<block_size; bi+=4) {
                    CONCAT(AVX_FUNC_NAME,_)(i00,j0,k0,bi,
                                            out, inL, inR, n);
                }
            }
        }
    }
}

#undef AVX_OP
#undef AVX_FUNC_NAME

