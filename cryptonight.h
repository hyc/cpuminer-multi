#ifndef __CRYPTONIGHT_H_INCLUDED
#define __CRYPTONIGHT_H_INCLUDED

#include <stddef.h>
#include "crypto/oaes_lib.h"
#include "miner.h"
#include "variant2_int_sqrt.h"

#define MEMORY         (1 << 21) /* 2 MiB */
#define ITER           (1 << 20)
#define AES_BLOCK_SIZE  16
#define AES_KEY_SIZE    32 /*16*/
#define INIT_SIZE_BLK   8
#define INIT_SIZE_BYTE (INIT_SIZE_BLK * AES_BLOCK_SIZE)	// 128

#pragma pack(push, 1)
union hash_state {
  uint8_t b[200];
  uint64_t w[25];
};
#pragma pack(pop)

#pragma pack(push, 1)
union cn_slow_hash_state {
    union hash_state hs;
    struct {
        uint8_t k[64];
        uint8_t init[INIT_SIZE_BYTE];
    };
};
#pragma pack(pop)

#ifdef USE_LOBOTOMIZED_AES

struct cryptonight_ctx {
    uint8_t long_state[MEMORY] __attribute((aligned(16)));
    union cn_slow_hash_state state;
    uint8_t text[INIT_SIZE_BYTE] __attribute((aligned(16)));
    uint8_t a[AES_BLOCK_SIZE] __attribute__((aligned(16)));
    uint8_t b[AES_BLOCK_SIZE * 2] __attribute__((aligned(16)));
    uint8_t c[AES_BLOCK_SIZE] __attribute__((aligned(16)));
    oaes_ctx* aes_ctx;
};

#else

struct cryptonight_ctx {
    uint8_t long_state[MEMORY] __attribute((aligned(16)));
    union cn_slow_hash_state state;
    uint8_t text[INIT_SIZE_BYTE] __attribute((aligned(16)));
    uint64_t a[AES_BLOCK_SIZE >> 3] __attribute__((aligned(16)));
    uint64_t b[AES_BLOCK_SIZE >> 2] __attribute__((aligned(16)));
    uint8_t c[AES_BLOCK_SIZE] __attribute__((aligned(16)));
    oaes_ctx* aes_ctx;
};

#endif

void do_blake_hash(const void* input, size_t len, char* output);
void do_groestl_hash(const void* input, size_t len, char* output);
void do_jh_hash(const void* input, size_t len, char* output);
void do_skein_hash(const void* input, size_t len, char* output);
void xor_blocks_dst(const uint8_t *restrict a, const uint8_t *restrict b, uint8_t *restrict dst);
int cryptonight_hash_ctx(void* output, const void* input, int inlen, struct cryptonight_ctx* ctx, int variant);
void keccak(const uint8_t *in, int inlen, uint8_t *md, int mdlen);
void keccakf(uint64_t st[25], int rounds);
extern void (* const extra_hashes[4])(const void *, size_t, char *);

#define VARIANT1_1(p) \
  do if (variant == 1) \
  { \
    const uint8_t tmp = ((const uint8_t*)(p))[11]; \
    static const uint32_t table = 0x75310; \
    const uint8_t index = (((tmp >> 3) & 6) | (tmp & 1)) << 1; \
    ((uint8_t*)(p))[11] = tmp ^ ((table >> index) & 0x30); \
  } while(0)

#define VARIANT1_INIT() \
  if (variant == 1 && inlen < 43) \
  { \
    return 0; \
  } \
  const uint64_t tweak1_2 = variant == 1 ? *((const uint64_t*) (((const uint8_t*)input) + 35)) ^ ctx->state.hs.w[24] : 0

#define VARIANT2_INIT64(b, state) \
  uint64_t division_result = 0; \
  uint64_t sqrt_result = 0; \
  do if (variant >= 2) \
  { \
    ((uint64_t*)b)[2] = state.hs.w[8] ^ state.hs.w[10]; \
    ((uint64_t*)b)[3] = state.hs.w[9] ^ state.hs.w[11]; \
    division_result = state.hs.w[12]; \
    sqrt_result = state.hs.w[13]; \
  } while (0)

#define VARIANT2_SHUFFLE_ADD_SSE2(base_ptr, offset) \
  do if (variant >= 2) \
  { \
    const __m128i chunk1 = _mm_load_si128((__m128i *)((base_ptr) + ((offset) ^ 0x10))); \
    const __m128i chunk2 = _mm_load_si128((__m128i *)((base_ptr) + ((offset) ^ 0x20))); \
    const __m128i chunk3 = _mm_load_si128((__m128i *)((base_ptr) + ((offset) ^ 0x30))); \
    _mm_store_si128((__m128i *)((base_ptr) + ((offset) ^ 0x10)), _mm_add_epi64(chunk3, b_x1)); \
    _mm_store_si128((__m128i *)((base_ptr) + ((offset) ^ 0x20)), _mm_add_epi64(chunk1, b_x)); \
    _mm_store_si128((__m128i *)((base_ptr) + ((offset) ^ 0x30)), _mm_add_epi64(chunk2, a_x)); \
  } while (0)

#define VARIANT2_PORTABLE_SHUFFLE_ADD(base_ptr, offset, a, b) \
  do if (variant >= 2) \
  { \
    uint64_t* chunk1 = (uint64_t*)((base_ptr) + ((offset) ^ 0x10)); \
    uint64_t* chunk2 = (uint64_t*)((base_ptr) + ((offset) ^ 0x20)); \
    uint64_t* chunk3 = (uint64_t*)((base_ptr) + ((offset) ^ 0x30)); \
    \
    const uint64_t chunk1_old[2] = { chunk1[0], chunk1[1] }; \
    \
    uint64_t b1[2]; \
    memcpy(b1, b + 16, 16); \
    chunk1[0] = chunk3[0] + b1[0]; \
    chunk1[1] = chunk3[1] + b1[1]; \
    \
    uint64_t a0[2]; \
    memcpy(a0, a, 16); \
    chunk3[0] = chunk2[0] + a0[0]; \
    chunk3[1] = chunk2[1] + a0[1]; \
    \
    uint64_t b0[2]; \
    memcpy(b0, b, 16); \
    chunk2[0] = chunk1_old[0] + b0[0]; \
    chunk2[1] = chunk1_old[1] + b0[1]; \
  } while (0)

#define VARIANT2_INTEGER_MATH_DIVISION_STEP(b, ptr) \
  ((uint64_t*)(b))[0] ^= division_result ^ (sqrt_result << 32); \
  { \
    const uint64_t dividend = ((uint64_t*)(ptr))[1]; \
    const uint32_t divisor = (((uint64_t*)(ptr))[0] + (uint32_t)(sqrt_result << 1)) | 0x80000001UL; \
    division_result = ((uint32_t)(dividend / divisor)) + \
                     (((uint64_t)(dividend % divisor)) << 32); \
  } \
  const uint64_t sqrt_input = ((uint64_t*)(ptr))[0] + division_result

#define VARIANT2_INTEGER_MATH_SSE2(b, ptr) \
  do if (variant >= 2) \
  { \
    VARIANT2_INTEGER_MATH_DIVISION_STEP(b, ptr); \
    VARIANT2_INTEGER_MATH_SQRT_STEP_SSE2(); \
    VARIANT2_INTEGER_MATH_SQRT_FIXUP(sqrt_result); \
  } while(0)

#define VARIANT2_2(hp_state, j) \
  do if (variant >= 2) \
  { \
    *(uint64_t*)(hp_state + (j ^ 0x10)) ^= hi; \
    *((uint64_t*)(hp_state + (j ^ 0x10)) + 1) ^= lo; \
    hi ^= *(uint64_t*)(hp_state + (j ^ 0x20)); \
    lo ^= *((uint64_t*)(hp_state + (j ^ 0x20)) + 1); \
  } while (0)

#endif
