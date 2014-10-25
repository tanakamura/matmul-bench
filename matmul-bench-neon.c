#if 0

#ifdef __ARM_NEON__

static NOINLINE void
neon_(unsigned int i00,
      unsigned int j0,
      unsigned int k0,
      unsigned int bi,
      float *__restrict out,
      const float *__restrict inL,
      const float *__restrict inR,
      unsigned int n, unsigned int pitch_f32)
{
    int i0 = i00+bi+0;

    float *outp0 = &out[i0*n+j0];

#define outp1 (outp0+n*1)

    float32x4_t *outp_0 = (float32x4_t*)outp0;
    float32x4_t *outp_1 = (float32x4_t*)outp1;

    __builtin_prefetch(outp1+n*1);
    __builtin_prefetch(outp1+n*2);

    float32x4_t vout0_0;
    float32x4_t vout0_1;
    float32x4_t vout0_2;
    float32x4_t vout0_3;

    float32x4_t vout1_0;
    float32x4_t vout1_1;
    float32x4_t vout1_2;
    float32x4_t vout1_3;


    if (k0==0) {
        vout0_0 = vdupq_n_f32(0);
        vout0_1 = vdupq_n_f32(0);
        vout0_2 = vdupq_n_f32(0);
        vout0_3 = vdupq_n_f32(0);
        vout1_0 = vdupq_n_f32(0);
        vout1_1 = vdupq_n_f32(0);
        vout1_2 = vdupq_n_f32(0);
        vout1_3 = vdupq_n_f32(0);
    } else {
        vout0_0 = outp_0[0];
        vout0_1 = outp_0[1];
        vout0_2 = outp_0[2];
        vout0_3 = outp_0[3];
        vout1_0 = outp_1[0];
        vout1_1 = outp_1[1];
        vout1_2 = outp_1[2];
        vout1_3 = outp_1[3];
    }

    const float *inRp1 = (float*)&inR[k0*pitch_f32+j0];
    const float *inRp1_pld = (float*)&inR[(k0+5)*pitch_f32+j0];

    const float *__restrict inL00_0 = (inL + (i0+0)*pitch_f32 + k0);
    const float *__restrict inL00_1 = (inL + (i0+1)*pitch_f32 + k0);

    __asm__ __volatile__ ("movs r4, #64\n"
                          ".p2align 4\n\t"
                          "1:\n\t"
                          "pld [%[inRp1_pld]]\n\t"
                          "vldmia %[inRp1], {q8-q11}\n\t"

                          "add %[inRp1], %[inRp1], %[pitch_f32]\n\t"
                          "vld1.32 {d24[], d25[]}, [%[inL00_0]]!\n\t"

                          "add %[inRp1_pld], %[inRp1_pld], %[pitch_f32]\n\t"
                          "vld1.32 {d26[], d27[]}, [%[inL00_1]]!\n\t"

                          "subs r4, r4, #1\n\t"

                          "vmla.f32 %q[vout0_0], q8, q12\n\t"
                          "vmla.f32 %q[vout0_1], q9, q12\n\t"
                          "vmla.f32 %q[vout0_2], q10, q12\n\t"
                          "vmla.f32 %q[vout0_3], q11, q12\n\t"
                          "vmla.f32 %q[vout1_0], q8, q13\n\t"
                          "vmla.f32 %q[vout1_1], q9, q13\n\t"
                          "vmla.f32 %q[vout1_2], q10, q13\n\t"
                          "vmla.f32 %q[vout1_3], q11, q13\n\t"
                          "bne 1b\n\t"
                          :[inL00_0]"+r"(inL00_0), [inL00_1]"+r"(inL00_1),
                           [vout0_0]"+w"(vout0_0),
                           [vout0_1]"+w"(vout0_1),
                           [vout0_2]"+w"(vout0_2),
                           [vout0_3]"+w"(vout0_3),
                           [vout1_0]"+w"(vout1_0),
                           [vout1_1]"+w"(vout1_1),
                           [vout1_2]"+w"(vout1_2),
                           [vout1_3]"+w"(vout1_3),
                           [inRp1]"+r"(inRp1),
                           [inRp1_pld]"+r"(inRp1_pld)
                          :[pitch_f32]"r"(pitch_f32*4)
                          :"r4", "q8", "q9", "q10", "q11", "q12", "q13");


    outp_0[0] = vout0_0;
    outp_0[1] = vout0_1;
    outp_0[2] = vout0_2;
    outp_0[3] = vout0_3;

    outp_1[0] = vout1_0;
    outp_1[1] = vout1_1;
    outp_1[2] = vout1_2;
    outp_1[3] = vout1_3;
}

static void
matmul_neon(float * __restrict out,
            const float* __restrict inL,
            const float* __restrict inR,
            unsigned int n,
            unsigned int pitch_f32)
{
    /* C=4x4x(2simd) register */
    unsigned int block_size_i = 32;
    unsigned int block_size_j = 16;
    unsigned int block_size_k = 64;
    int i00;

#pragma omp parallel for schedule(dynamic)
    for (i00=0; i00<n; i00+=block_size_i) {
        for (int j0=0; j0<n; j0+=block_size_j) {
            for (int k0=0; k0<n; k0+=block_size_k) {
                for (int bi=0; bi<block_size_i; bi+=2) {
                    neon_(i00, j0, k0, bi, out, inL, inR, n, pitch_f32);
                }
            }
        }
    }
}
#endif



#endif