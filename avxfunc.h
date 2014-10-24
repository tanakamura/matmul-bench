#ifndef _mm256_set_m128
#define _mm256_set_m128(hi,lo) \
    _mm256_insertf128_ps(_mm256_castps128_ps256(lo), (hi), 0x1)
#endif


static inline void
CONCAT(AVX_FUNC_NAME,_)(unsigned long i00,
                        unsigned long j0,
                        unsigned long k0,
                        unsigned long bi,
                        float * __restrict out,
                        const float* __restrict inL,
                        const float* __restrict inR,
                        unsigned long n,
                        unsigned long pitch_f32)
{
    long i0 = i00+bi+0;
    float *outp0 = &out[i0*n+j0];

#define outp1 (outp0+n*1)
#define outp2 (outp0+n*2)

    _mm_prefetch((const char*)(outp0), _MM_HINT_T0);
    _mm_prefetch((const char*)(outp1), _MM_HINT_T0);
    _mm_prefetch((const char*)(outp2), _MM_HINT_T0);

    __m256 *outp8_0 = (__m256*)outp0;
    __m256 *outp8_1 = (__m256*)outp1;
    __m256 *outp8_2 = (__m256*)outp2;

    __m256 vout0_0 = _mm256_setzero_ps();
    __m256 vout0_1 = _mm256_setzero_ps();
    __m256 vout0_2 = _mm256_setzero_ps();
    __m256 vout0_3 = _mm256_setzero_ps();

    __m256 vout1_0 = _mm256_setzero_ps();
    __m256 vout1_1 = _mm256_setzero_ps();
    __m256 vout1_2 = _mm256_setzero_ps();
    __m256 vout1_3 = _mm256_setzero_ps();

    __m256 vout2_0 = _mm256_setzero_ps();
    __m256 vout2_1 = _mm256_setzero_ps();
    __m256 vout2_2 = _mm256_setzero_ps();
    __m256 vout2_3 = _mm256_setzero_ps();

#define AVX_K8(K)                                               \
    {                                                           \
        const float *inL00 = inL0;                              \
                                                                \
        _mm_prefetch((const char*)(inR0 + pitch_f32*4), _MM_HINT_T0);   \
                                                                \
        lik0_8 = _mm256_set1_ps(*inL00);                         \
        inL00+=pitch_f32;                                        \
        lik1_8 = _mm256_set1_ps(*inL00);                        \
        inL00+=pitch_f32;                                        \
        lik2_8 = _mm256_set1_ps(*inL00);                         \
                                                                \
        vr0 = _mm256_load_ps(&inR0[0]);                         \
        AVX_OP(0,0);                                            \
        AVX_OP(1,0);                                            \
        AVX_OP(2,0);                                            \
                                                                \
        vr1 = _mm256_load_ps(&inR0[8]);                         \
        AVX_OP(0,1);                                            \
        AVX_OP(1,1);                                            \
        AVX_OP(2,1);                                            \
                                                                \
        vr2 = _mm256_load_ps(&inR0[16]);                        \
        AVX_OP(0,2);                                            \
        AVX_OP(1,2);                                            \
        AVX_OP(2,2);                                            \
                                                                \
        vr3 = _mm256_load_ps(&inR0[24]);                        \
        AVX_OP(0,3);                                            \
        AVX_OP(1,3);                                            \
        AVX_OP(2,3);                                            \
        inL0++;                                                 \
        inR0 += pitch_f32;                                      \
    }


    __m256 lik0_8, lik1_8, lik2_8;
    __m256 vr0, vr1, vr2, vr3;

#define LC(N) AVX_K8(N);

    const float *inL0 = &inL[i0*pitch_f32+k0];
    const float *inR0 = &inR[k0*pitch_f32+j0];

    for (int i=0; i<8; i++) {
        M16(LC);
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

        outp8_2[0] = vout2_0;
        outp8_2[1] = vout2_1;
        outp8_2[2] = vout2_2;
        outp8_2[3] = vout2_3;
    } else {
        outp8_0[0] = _mm256_add_ps(outp8_0[0], vout0_0);
        outp8_0[1] = _mm256_add_ps(outp8_0[1], vout0_1);
        outp8_0[2] = _mm256_add_ps(outp8_0[2], vout0_2);
        outp8_0[3] = _mm256_add_ps(outp8_0[3], vout0_3);

        outp8_1[0] = _mm256_add_ps(outp8_1[0], vout1_0);
        outp8_1[1] = _mm256_add_ps(outp8_1[1], vout1_1);
        outp8_1[2] = _mm256_add_ps(outp8_1[2], vout1_2);
        outp8_1[3] = _mm256_add_ps(outp8_1[3], vout1_3);

        outp8_2[0] = _mm256_add_ps(outp8_2[0], vout2_0);
        outp8_2[1] = _mm256_add_ps(outp8_2[1], vout2_1);
        outp8_2[2] = _mm256_add_ps(outp8_2[2], vout2_2);
        outp8_2[3] = _mm256_add_ps(outp8_2[3], vout2_3);
    }

}


static void
AVX_FUNC_NAME(float * __restrict out,
              const float* __restrict inL,
              const float* __restrict inR,
              unsigned long n,
              unsigned long pitch_f32)
{
    unsigned long block_size_i = 48;
    unsigned long block_size_j = 32;
    unsigned long block_size_k = 128;

    long i00;

#pragma omp parallel for schedule(static)
    for (i00=0; i00<n; i00+=block_size_i) {
        for (long j0=0; j0<n; j0+=block_size_j) {
            for (long k0=0; k0<n; k0+=block_size_k) {
                for (long bi=0; bi<block_size_i; bi+=3) {
                    CONCAT(AVX_FUNC_NAME,_)(i00,j0,k0,bi,
                                            out, inL, inR, n, pitch_f32);
                }
            }
        }
    }
}

#undef AVX_OP
#undef AVX_FUNC_NAME

