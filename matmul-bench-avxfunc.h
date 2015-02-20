#ifndef _mm256_set_m128
#define _mm256_set_m128(hi,lo) \
    _mm256_insertf128_ps(_mm256_castps128_ps256(lo), (hi), 0x1)
#endif

#ifdef MAT_4x2
#define AVX_OP_2(J) 
#else
#define AVX_OP_2(J) AVX_OP(2,J)
#endif

static NOINLINE W32_ALIGN_ARG_POINTER void
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

#ifdef MAT_4x3
#define outp2 (outp0+n*2)
#endif

    _mm_prefetch((const char*)(outp0), _MM_HINT_T0);
    _mm_prefetch((const char*)(outp1), _MM_HINT_T0);
#ifdef MAT_4x3
    _mm_prefetch((const char*)(outp2), _MM_HINT_T0);
#endif

    __m256 *outp8_0 = (__m256*)outp0;
    __m256 *outp8_1 = (__m256*)outp1;

#ifdef MAT_4x3
    __m256 *outp8_2 = (__m256*)outp2;
#endif

    __m256 vout0_0 = _mm256_setzero_ps();
    __m256 vout0_1 = _mm256_setzero_ps();
    __m256 vout0_2 = _mm256_setzero_ps();
    __m256 vout0_3 = _mm256_setzero_ps();

    __m256 vout1_0 = _mm256_setzero_ps();
    __m256 vout1_1 = _mm256_setzero_ps();
    __m256 vout1_2 = _mm256_setzero_ps();
    __m256 vout1_3 = _mm256_setzero_ps();

#ifdef MAT_4x3
    __m256 vout2_0 = _mm256_setzero_ps();
    __m256 vout2_1 = _mm256_setzero_ps();
    __m256 vout2_2 = _mm256_setzero_ps();
    __m256 vout2_3 = _mm256_setzero_ps();
#endif

#ifdef MAT_4x3
#define LOAD_LIK2_8() lik2_8 = _mm256_set1_ps(*inL00);
#else
#define LOAD_LIK2_8() 
#endif

#define AVX_K8(K)                                               \
    {                                                           \
        const float *inL00 = inL0;                              \
                                                                \
        _mm_prefetch((const char*)(inR0 + pitch_f32*4), _MM_HINT_T0);   \
                                                                \
        lik0_8 = _mm256_set1_ps(*inL00);                         \
        inL00+=pitch_f32;                                        \
        vr0 = _mm256_load_ps(&inR0[0]);                         \
                                                                \
        AVX_OP(0,0);                                            \
        lik1_8 = _mm256_set1_ps(*inL00);                         \
        inL00+=pitch_f32;                                        \
        AVX_OP(1,0);                                            \
        LOAD_LIK2_8();                                           \
        AVX_OP_2(0);                                             \
                                                                \
        vr1 = _mm256_load_ps(&inR0[8]);                         \
        AVX_OP(0,1);                                            \
        AVX_OP(1,1);                                            \
        AVX_OP_2(1);                                            \
                                                                \
        vr2 = _mm256_load_ps(&inR0[16]);                        \
        AVX_OP(0,2);                                            \
        AVX_OP(1,2);                                            \
        AVX_OP_2(2);                                            \
                                                                \
        vr3 = _mm256_load_ps(&inR0[24]);                        \
        AVX_OP(0,3);                                            \
        AVX_OP(1,3);                                            \
        AVX_OP_2(3);                                            \
        inL0++;                                                 \
        inR0 += pitch_f32;                                      \
    }


    __m256 lik0_8, lik1_8;
#ifdef MAT_4x3
    __m256 lik2_8;
#endif
    __m256 vr0, vr1, vr2, vr3;

#define LC(N) AVX_K8(N);

    const float *inL0 = &inL[i0*pitch_f32+k0];
    const float *inR0 = &inR[k0*pitch_f32+j0];

    //_mm_prefetch(inR0 + pitch_f32*1, _MM_HINT_T0);
    //_mm_prefetch(inR0 + pitch_f32*2, _MM_HINT_T0);
    //_mm_prefetch(inR0 + pitch_f32*3, _MM_HINT_T0);

    for (int i=0; i<32; i++) {
        _mm_prefetch(inL0 + 16, _MM_HINT_T0);
        M4(LC);
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
#ifdef MAT_4x3
        outp8_2[0] = vout2_0;
        outp8_2[1] = vout2_1;
        outp8_2[2] = vout2_2;
        outp8_2[3] = vout2_3;
#endif
    } else {
        outp8_0[0] = _mm256_add_ps(outp8_0[0], vout0_0);
        outp8_0[1] = _mm256_add_ps(outp8_0[1], vout0_1);
        outp8_0[2] = _mm256_add_ps(outp8_0[2], vout0_2);
        outp8_0[3] = _mm256_add_ps(outp8_0[3], vout0_3);

        outp8_1[0] = _mm256_add_ps(outp8_1[0], vout1_0);
        outp8_1[1] = _mm256_add_ps(outp8_1[1], vout1_1);
        outp8_1[2] = _mm256_add_ps(outp8_1[2], vout1_2);
        outp8_1[3] = _mm256_add_ps(outp8_1[3], vout1_3);
#ifdef MAT_4x3
        outp8_2[0] = _mm256_add_ps(outp8_2[0], vout2_0);
        outp8_2[1] = _mm256_add_ps(outp8_2[1], vout2_1);
        outp8_2[2] = _mm256_add_ps(outp8_2[2], vout2_2);
        outp8_2[3] = _mm256_add_ps(outp8_2[3], vout2_3);
#endif
    }

}

static void
CONCAT(AVX_FUNC_NAME,thread)(struct MatmulBenchParam *p,
                             unsigned long i_start,
                             unsigned long i_end,
                             unsigned int thread_id)
{
    float * __restrict out = p->out;
    unsigned long n = p->n;
    const float * __restrict inL_plus1line = p->inL_plus1line;
    const float * __restrict inR_plus1line = p->inR_plus1line;

    unsigned long pitch_byte = p->pitch_byte;

#ifdef MAT_4x3
    unsigned long block_size_i = 48;
#else
    unsigned long block_size_i = 64;
#endif
    unsigned long block_size_j = 32;
    unsigned long block_size_k = 128;

    for (unsigned long i00=i_start; i00<i_end; i00+=block_size_i) {
        for (long j0=0; j0<n; j0+=block_size_j) {
            for (long k0=0; k0<n; k0+=block_size_k) {
#ifdef MAT_4x3
                unsigned int i_inc = 3;
#else
                unsigned int i_inc = 2;

#endif

                for (long bi=0; bi<block_size_i; bi+=i_inc) {
                    CONCAT(AVX_FUNC_NAME,_)(i00,j0,k0,bi,
                                            out, inL_plus1line, inR_plus1line, n, pitch_byte/4);
                }
            }
        }
    }
}



static void
AVX_FUNC_NAME(struct MatmulBenchParam *p)
{
#ifdef MAT_4x3
    unsigned long block_size_i = 48;
#else
    unsigned long block_size_i = 64;
#endif

    matmul_bench_thread_call(p, p->i_block_size*block_size_i, p->n, CONCAT(AVX_FUNC_NAME,thread));
}

#undef AVX_OP
#undef AVX_FUNC_NAME

