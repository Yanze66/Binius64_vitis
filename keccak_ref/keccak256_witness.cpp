#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define N_CONST 25
#define N_WITNESS 5
#define N_INTERNAL 721
#define COMMITTED 1024
#define N_SCRATCH 2057
#define TOTAL (COMMITTED + N_SCRATCH)

typedef struct {
    uint32_t idx;
} Wire;

uint64_t value_vec[TOTAL];

// ---------- builder state ----------
uint32_t next_const = 0;
uint32_t next_witness = N_CONST;
uint32_t next_internal = N_CONST + N_WITNESS;

// ---------- constant pool ----------
typedef struct {
    uint64_t value;
    uint32_t idx;
} ConstEntry;

ConstEntry const_pool[128];
size_t const_pool_len = 0;

// ---------- helpers ----------
Wire add_constant(uint64_t v) {
    for (size_t i = 0; i < const_pool_len; i++) {
        if (const_pool[i].value == v) {
            return (Wire){ const_pool[i].idx };
        }
    }
    uint32_t idx = next_const++;
    value_vec[idx] = v;
    const_pool[const_pool_len++] = (ConstEntry){ v, idx };
    return (Wire){ idx };
}

Wire add_witness(void) {
    return (Wire){ next_witness++ };
}

Wire add_internal(uint64_t v) {
    uint32_t idx = next_internal++;
    value_vec[idx] = v;
    return (Wire){ idx };
}

// ---------- gates ----------
Wire bxor(Wire a, Wire b) {
    return add_internal(value_vec[a.idx] ^ value_vec[b.idx]);
}

Wire band(Wire a, Wire b) {
    return add_internal(value_vec[a.idx] & value_vec[b.idx]);
}

Wire bnot(Wire a) {
    return add_internal(~value_vec[a.idx]);
}

static inline uint64_t rotl64(uint64_t x, uint32_t r) {
    return (x << r) | (x >> (64 - r));
}

Wire rotl(Wire a, uint32_t r) {
    return add_internal(rotl64(value_vec[a.idx], r));
}

Wire fax(Wire na, Wire b, Wire c) {
    // (~a & b) ^ c
    return add_internal((value_vec[na.idx] & value_vec[b.idx]) ^ value_vec[c.idx]);
}

Wire bxor_multi(Wire *w, size_t n) {
    uint64_t acc = 0;
    for (size_t i = 0; i < n; i++) {
        acc ^= value_vec[w[i].idx];
    }
    return add_internal(acc);
}

// ---------- keccak constants ----------
static const uint64_t RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL,
    0x800000000000808AULL, 0x8000000080008000ULL,
    0x000000000000808BULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008AULL, 0x0000000000000088ULL,
    0x0000000080008009ULL, 0x000000008000000AULL,
    0x000000008000808BULL, 0x800000000000008BULL,
    0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800AULL, 0x800000008000000AULL,
    0x8000000080008081ULL, 0x8000000000008080ULL,
    0x0000000080000001ULL, 0x8000000080008008ULL
};

static const uint32_t R[25] = {
     0,  1, 62, 28, 27,
    36, 44,  6, 55, 20,
     3, 10, 43, 25, 39,
    41, 45, 15, 21,  8,
    18,  2, 61, 56, 14
};

static inline uint32_t idx(uint32_t x, uint32_t y) {
    return x + 5*y;
}

// ---------- permutation ----------
void keccak_f1600(Wire state[25]) {
    for (int round = 0; round < 24; round++) {
        // theta
        Wire C[5];
        for (int x = 0; x < 5; x++) {
            Wire col[5];
            for (int y = 0; y < 5; y++) col[y] = state[idx(x,y)];
            C[x] = bxor_multi(col, 5);
        }

        Wire D[5];
        for (int x = 0; x < 5; x++) {
            D[x] = bxor(C[(x+4)%5], rotl(C[(x+1)%5], 1));
        }

        for (int y = 0; y < 5; y++) {
            for (int x = 0; x < 5; x++) {
                state[idx(x,y)] = bxor(state[idx(x,y)], D[x]);
            }
        }

        // rho + pi
        Wire tmp[25];
        memcpy(tmp, state, sizeof(tmp));
        for (int y = 0; y < 5; y++) {
            for (int x = 0; x < 5; x++) {
                uint32_t r = R[idx(x,y)];
                if (r != 0) {
                    tmp[idx(y,(2*x+3*y)%5)] = rotl(state[idx(x,y)], r);
                }
            }
        }
        memcpy(state, tmp, sizeof(tmp));

        // chi
        for (int y = 0; y < 5; y++) {
            Wire a0 = state[idx(0,y)];
            Wire a1 = state[idx(1,y)];
            Wire a2 = state[idx(2,y)];
            Wire a3 = state[idx(3,y)];
            Wire a4 = state[idx(4,y)];

            state[idx(0,y)] = fax(bnot(a1), a2, a0);
            state[idx(1,y)] = fax(bnot(a2), a3, a1);
            state[idx(2,y)] = fax(bnot(a3), a4, a2);
            state[idx(3,y)] = fax(bnot(a4), a0, a3);
            state[idx(4,y)] = fax(bnot(a0), a1, a4);
        }

        // iota
        Wire rc = add_constant(RC[round]);
        state[0] = bxor(state[0], rc);
    }
}


// ====== 来自你给的测试 ======
constexpr size_t LEN = 8;
uint8_t msg[LEN] = { 178, 96, 184, 161, 3, 67, 191, 90 };

// golden digest (Rust sha3 crate)
uint8_t expected_digest_bytes[32] = {
    0xf3,0xa4,0x33,0x17,0x36,0xe3,0x9b,0x61,
    0xb5,0xea,0x9c,0x4f,0x07,0x66,0x12,0xe3,
    0xba,0xb7,0xd4,0xcf,0x36,0xce,0x5a,0xf1,
    0x86,0x31,0xd7,0xad,0xc5,0xde,0x69,0xe3
};

// ====== helper ======
uint64_t le_bytes_to_u64(const uint8_t *b, size_t n) {
    uint64_t v = 0;
    for (size_t i = 0; i < n; i++) {
        v |= (uint64_t)b[i] << (8 * i);
    }
    return v;
}

int main() {
    printf("=== Keccak256 CIRCUIT test (LEN = 8) ===\n");

    // --------------------------------------------------
    // 1. 创建 message witness wires
    // --------------------------------------------------
    Wire message_wires[1]; // LEN=8 → 1 word
    message_wires[0] = add_witness();

    // --------------------------------------------------
    // 2. 创建 expected digest witness wires
    // --------------------------------------------------
    Wire expected_digest_wires[4];
    for (int i = 0; i < 4; i++) {
        expected_digest_wires[i] = add_witness();
    }

    // --------------------------------------------------
    // 3. padding（len_bytes.is_multiple_of(8) 分支）
    // --------------------------------------------------
    // padded_message = [message_word, 0x01, 0, 0, ..., last ^= 0x80<<56]
    Wire padded[17]; // 1 block = 17 words

    padded[0] = message_wires[0];
    padded[1] = add_constant(0x01ULL);

    Wire zero = add_constant(0);
    for (int i = 2; i < 17; i++) {
        padded[i] = zero;
    }

    // XOR 0x80 into last word MSB
    Wire last_mask = add_constant(0x80ULL << 56);
    padded[16] = bxor(padded[16], last_mask);

    // --------------------------------------------------
    // 4. 初始化 state = 0
    // --------------------------------------------------
    Wire state[25];
    for (int i = 0; i < 25; i++) {
        state[i] = zero;
    }

    // --------------------------------------------------
    // 5. absorb
    // --------------------------------------------------
    for (int i = 0; i < 17; i++) {
        state[i] = bxor(state[i], padded[i]);
    }

    // --------------------------------------------------
    // 6. keccak-f[1600]
    // --------------------------------------------------
    keccak_f1600(state);

    // --------------------------------------------------
    // 7. computed digest = state[0..3]
    // --------------------------------------------------
    Wire computed_digest[4] = {
        state[0], state[1], state[2], state[3]
    };

    // --------------------------------------------------
    // 8. “assert_eq” → 在 C 里等价为 runtime check
    // --------------------------------------------------
    for (int i = 0; i < 4; i++) {
        uint64_t cd = value_vec[computed_digest[i].idx];
        uint64_t ed = le_bytes_to_u64(expected_digest_bytes + 8*i, 8);
        if (cd != ed) {
            printf("ASSERT FAIL digest[%d]\n", i);
            return 1;
        }
    }

    printf("Digest matches golden ✔\n");

    // --------------------------------------------------
    // 9. 填 witness（等价于 WitnessFiller）
    // --------------------------------------------------
    value_vec[message_wires[0].idx] =
        le_bytes_to_u64(msg, 8);

    for (int i = 0; i < 4; i++) {
        value_vec[expected_digest_wires[i].idx] =
            le_bytes_to_u64(expected_digest_bytes + 8*i, 8);
    }

    // --------------------------------------------------
    // 10. dump 完整 witness
    // --------------------------------------------------
    printf("\n=== FULL WITNESS DUMP ===\n");
    for (int i = 0; i < TOTAL; i++) {
        printf("W[%04d] = 0x%016llx\n", i,
               (unsigned long long)value_vec[i]);
    }

    return 0;
}
