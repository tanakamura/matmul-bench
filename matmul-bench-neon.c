#include "matmul-bench-common.h"

#define BLOCK_SIZE_I 64
#define BLOCK_SIZE_J 16
#define BLOCK_SIZE_K 128

static inline void
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

    vout0_0 = vdupq_n_f32(0);
    vout0_1 = vdupq_n_f32(0);
    vout0_2 = vdupq_n_f32(0);
    vout0_3 = vdupq_n_f32(0);
    vout1_0 = vdupq_n_f32(0);
    vout1_1 = vdupq_n_f32(0);
    vout1_2 = vdupq_n_f32(0);
    vout1_3 = vdupq_n_f32(0);

    const float *inRp1 = (float*)&inR[k0*pitch_f32+j0];

    //const float *__restrict inL_orig = inL;
    const float *__restrict inR_orig = inRp1;

    const float *__restrict inL00_0 = (inL + (i0+0)*pitch_f32 + k0);
    const float *__restrict inL00_1 = (inL + (i0+1)*pitch_f32 + k0);
    const float *__restrict inL_next = inL00_1 + pitch_f32;

#define K_LOOP_BODY(PLD_L0,PLD_L1,PLD_R0,PLD_R1,PLD_R2,PLD_R3,CMP)      \
    ""                                                                  \
        "vld1.32 {d0,d1}, [%[inL00_0]:64]!\n\t"                         \
        "vld1.32 {d2,d3}, [%[inL00_1]:64]!\n\t"                         \
        PLD_R0                                                          \
        PLD_L0                                                          \
        "vldmia %[inRp1], {q8-q11}\n\t"                                 \
                                                                        \
        "add %[inRp1], %[inRp1], %[pitch_f32]\n\t"                      \
                                                                        \
        "vmla.f32 %q[vout0_0], q8, d0[0]\n\t"                           \
        "vmla.f32 %q[vout1_0], q8, d2[0]\n\t"                           \
        "vmla.f32 %q[vout0_1], q9, d0[0]\n\t"                           \
        "vmla.f32 %q[vout1_1], q9, d2[0]\n\t"                           \
        "vmla.f32 %q[vout0_2], q10, d0[0]\n\t"                          \
        "vmla.f32 %q[vout1_2], q10, d2[0]\n\t"                          \
        "vmla.f32 %q[vout0_3], q11, d0[0]\n\t"                          \
        "vmla.f32 %q[vout1_3], q11, d2[0]\n\t"                          \
                                                                        \
        PLD_R1                                                          \
        "vldmia %[inRp1], {q8-q11}\n\t"                                 \
        "add %[inRp1], %[inRp1], %[pitch_f32]\n\t"                      \
                                                                        \
        "vmla.f32 %q[vout0_0], q8, d0[1]\n\t"                           \
        "vmla.f32 %q[vout1_0], q8, d2[1]\n\t"                           \
        "vmla.f32 %q[vout0_1], q9, d0[1]\n\t"                           \
        "vmla.f32 %q[vout1_1], q9, d2[1]\n\t"                           \
        "vmla.f32 %q[vout0_2], q10, d0[1]\n\t"                          \
        "vmla.f32 %q[vout1_2], q10, d2[1]\n\t"                          \
        "vmla.f32 %q[vout0_3], q11, d0[1]\n\t"                          \
        "vmla.f32 %q[vout1_3], q11, d2[1]\n\t"                          \
                                                                        \
        PLD_R2                                                          \
        PLD_L1                                                          \
        "vldmia %[inRp1], {q8-q11}\n\t"                                 \
        "add %[inRp1], %[inRp1], %[pitch_f32]\n\t"                      \
                                                                        \
        "vmla.f32 %q[vout0_0], q8, d1[0]\n\t"                           \
        "vmla.f32 %q[vout1_0], q8, d3[0]\n\t"                           \
        "vmla.f32 %q[vout0_1], q9, d1[0]\n\t"                           \
        "vmla.f32 %q[vout1_1], q9, d3[0]\n\t"                           \
        "vmla.f32 %q[vout0_2], q10, d1[0]\n\t"                          \
        "vmla.f32 %q[vout1_2], q10, d3[0]\n\t"                          \
        "vmla.f32 %q[vout0_3], q11, d1[0]\n\t"                          \
        "vmla.f32 %q[vout1_3], q11, d3[0]\n\t"                          \
                                                                        \
        PLD_R3                                                          \
        "vldmia %[inRp1], {q8-q11}\n\t"                                 \
        "add %[inRp1], %[inRp1], %[pitch_f32]\n\t"                      \
                                                                        \
        "vmla.f32 %q[vout0_0], q8, d1[1]\n\t"                           \
        "vmla.f32 %q[vout1_0], q8, d3[1]\n\t"                           \
        "vmla.f32 %q[vout0_1], q9, d1[1]\n\t"                           \
        "vmla.f32 %q[vout1_1], q9, d3[1]\n\t"                           \
        CMP                                                             \
        "vmla.f32 %q[vout0_2], q10, d1[1]\n\t"                          \
        "vmla.f32 %q[vout1_2], q10, d3[1]\n\t"                          \
        "vmla.f32 %q[vout0_3], q11, d1[1]\n\t"                          \
        "vmla.f32 %q[vout1_3], q11, d3[1]\n\t"

    __asm__ __volatile__ (".p2align 3\n\t"
                          "1:\n\t"
                          K_LOOP_BODY("pld [%[inL00_0], #64]\n\t",
                                      "pld [%[inL00_1], #64]\n\t",
                                      "pld [%[inRp1], %[pld_offset]]\n\t",
                                      "pld [%[inRp1], %[pld_offset]]\n\t",
                                      "pld [%[inRp1], %[pld_offset]]\n\t",
                                      "pld [%[inRp1], %[pld_offset]]\n\t",
                                      "cmp %[inRp1], %[inRp_end]\n\t")

                          "bne 1b\n\t"

                          K_LOOP_BODY("pld [%[inL00_0], #64]\n\t",
                                      "pld [%[inL00_1], #64]\n\t",
                                      "pld [%[inR_orig]]\n\t",
                                      "pld [%[inR_orig], %[pitch_f32]]\n\t",
                                      "pld [%[inR_orig], %[pitch_f32], lsl # 1]\n\t",
                                      "\n\t",
                                      "\n\t")

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
                           [inR_orig]"+r"(inR_orig),
                           [inL_next]"+r"(inL_next)
                          :[pitch_f32]"r"(pitch_f32*4), [inRp_end]"r"(inRp1 + (BLOCK_SIZE_K-4)*pitch_f32), [pld_offset]"r"(pitch_f32*4*4)
                          :"d0", "d1", "d2", "d3", "q8", "q9", "q10", "q11");

#undef K_LOOP_BODY

    if (k0 == 0) {
        outp_0[0] = vout0_0;
        outp_0[1] = vout0_1;
        outp_0[2] = vout0_2;
        outp_0[3] = vout0_3;

        outp_1[0] = vout1_0;
        outp_1[1] = vout1_1;
        outp_1[2] = vout1_2;
        outp_1[3] = vout1_3;
    } else {
        outp_0[0] = vaddq_f32(outp_0[0], vout0_0);
        outp_0[1] = vaddq_f32(outp_0[1], vout0_1);
        outp_0[2] = vaddq_f32(outp_0[2], vout0_2);
        outp_0[3] = vaddq_f32(outp_0[3], vout0_3);

        outp_1[0] = vaddq_f32(outp_1[0], vout1_0);
        outp_1[1] = vaddq_f32(outp_1[1], vout1_1);
        outp_1[2] = vaddq_f32(outp_1[2], vout1_2);
        outp_1[3] = vaddq_f32(outp_1[3], vout1_3);
    }
}

static void
neon_thread(struct MatmulBenchParam *p,
            unsigned long i_start,
            unsigned long i_end,
            unsigned int thread_id)
{
    float * __restrict out = p->out;
    const float * __restrict inL_plus1line = p->inL_plus1line;
    const float * __restrict inR_plus1line = p->inR_plus1line;
    unsigned int n = p->n;
    unsigned int pitch_byte = p->pitch_byte;

    unsigned int block_size_i = 32;
    unsigned int block_size_j = 16;
    unsigned int block_size_k = 128;

    for (unsigned long i00=i_start; i00<i_end; i00+=block_size_i) {
        for (int j0=0; j0<n; j0+=block_size_j) {
            for (int k0=0; k0<n; k0+=block_size_k) {
                for (int bi=0; bi<block_size_i; bi+=2) {
                    neon_(i00, j0, k0, bi, out, inL_plus1line, inR_plus1line, n, pitch_byte/4);
                }
            }
        }
    }
}

static void
neon_run(struct MatmulBenchParam *p)
{
    unsigned int block_size_i = 32;

    matmul_bench_thread_call(p, p->i_block_size*block_size_i, p->n, neon_thread);

}

static const struct MatmulBenchTest neon = MATMULBENCH_TEST_INITIALIZER("neon", neon_run, 128);

void
matmulbench_init_neon(struct MatmulBench *b, struct npr_varray *test_set)
{
    VA_PUSH(struct MatmulBenchTest, test_set, neon);
    //VA_PUSH(struct MatmulBenchTest, test_set, neon_4x3);
}
