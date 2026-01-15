// keccak_witness_aligned.c
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* =========================================================
 *  Layout (MUST match Rust dump)
 * ========================================================= */
#define N_CONST     25
#define N_WITNESS   5
#define COMMITTED   1024   // only committed, ignore scratch

typedef struct { uint32_t idx; } Wire;
static uint64_t value_vec[COMMITTED];
static uint32_t next_internal = N_CONST + N_WITNESS;

/* =========================================================
 *  Fixed constant table (EXACT order from Rust dump)
 * ========================================================= */
static const uint64_t CONST_TABLE[N_CONST] = {
    0x0000000000000000ULL,
    0x0000000000000001ULL,
    0x0000000000000088ULL,
    0x000000000000008AULL,
    0x000000000000800AULL,
    0x0000000000008082ULL,
    0x000000000000808BULL,
    0x0000000080000001ULL,
    0x000000008000000AULL,
    0x0000000080008009ULL,
    0x000000008000808BULL,
    0x8000000000000000ULL,
    0x8000000000000080ULL,
    0x800000000000008BULL,
    0x8000000000008002ULL,
    0x8000000000008003ULL,
    0x8000000000008009ULL,
    0x8000000000008080ULL,
    0x8000000000008089ULL,
    0x800000000000808AULL,
    0x800000008000000AULL,
    0x8000000080008000ULL,
    0x8000000080008008ULL,
    0x8000000080008081ULL,
    0xFFFFFFFFFFFFFFFFULL
};

/* =========================================================
 *  Builder helpers
 * ========================================================= */
static inline Wire constant_wire(uint32_t i) {
    assert(i < N_CONST);
    return (Wire){ i };
}

static inline Wire add_internal(uint64_t v) {
    assert(next_internal < COMMITTED);
    Wire w = { next_internal++ };
    value_vec[w.idx] = v;
    return w;
}

/* =========================================================
 *  Gates (1 gate = 1 internal)
 * ========================================================= */
static inline Wire bxor(Wire a, Wire b) {
    return add_internal(value_vec[a.idx] ^ value_vec[b.idx]);
}

static inline Wire bnot(Wire a) {
    return add_internal(~value_vec[a.idx]);
}

static inline uint64_t rotl64(uint64_t x, uint32_t r) {
    return (x << r) | (x >> (64 - r));
}

static inline Wire rotl(Wire a, uint32_t r) {
    return add_internal(rotl64(value_vec[a.idx], r));
}

static inline Wire fax(Wire na, Wire b, Wire c) {
    return add_internal((value_vec[na.idx] & value_vec[b.idx]) ^ value_vec[c.idx]);
}

static inline Wire bxor_multi(Wire *w, int n) {
    uint64_t acc = 0;
    for (int i = 0; i < n; i++) acc ^= value_vec[w[i].idx];
    return add_internal(acc);
}

/* =========================================================
 *  Keccak constants
 * ========================================================= */
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

static inline uint32_t IDX(uint32_t x, uint32_t y) {
    return x + 5*y;
}

/* =========================================================
 *  Keccak-f[1600] (witness-aligned)
 * ========================================================= */
void keccak_f1600(Wire state[25]) {
    for (int round = 0; round < 24; round++) {

        // theta
        Wire C[5];
        for (int x = 0; x < 5; x++) {
            Wire col[5];
            for (int y = 0; y < 5; y++) col[y] = state[IDX(x,y)];
            C[x] = bxor_multi(col, 5);
        }

        Wire D[5];
        for (int x = 0; x < 5; x++) {
            D[x] = bxor(C[(x+4)%5], rotl(C[(x+1)%5], 1));
        }

        for (int y = 0; y < 5; y++)
            for (int x = 0; x < 5; x++)
                state[IDX(x,y)] = bxor(state[IDX(x,y)], D[x]);

        // rho + pi
        Wire tmp[25];
        for (int i = 0; i < 25; i++) tmp[i] = state[0];

        for (int y = 0; y < 5; y++)
            for (int x = 0; x < 5; x++) {
                uint32_t r = R[IDX(x,y)];
                if (r)
                    tmp[IDX(y,(2*x+3*y)%5)] = rotl(state[IDX(x,y)], r);
            }

        for (int i = 0; i < 25; i++) state[i] = tmp[i];

        // chi
        for (int y = 0; y < 5; y++) {
            Wire a0 = state[IDX(0,y)];
            Wire a1 = state[IDX(1,y)];
            Wire a2 = state[IDX(2,y)];
            Wire a3 = state[IDX(3,y)];
            Wire a4 = state[IDX(4,y)];

            state[IDX(0,y)] = fax(bnot(a1), a2, a0);
            state[IDX(1,y)] = fax(bnot(a2), a3, a1);
            state[IDX(2,y)] = fax(bnot(a3), a4, a2);
            state[IDX(3,y)] = fax(bnot(a4), a0, a3);
            state[IDX(4,y)] = fax(bnot(a0), a1, a4);
        }

        // iota
        state[0] = bxor(state[0], constant_wire(1 + round));
    }

    // final internal (matches Rust)
    add_internal(0);
}

/* =========================================================
 *  Utils
 * ========================================================= */
static uint64_t le_bytes_to_u64(const uint8_t *b, int n) {
    uint64_t v = 0;
    for (int i = 0; i < n; i++)
        v |= (uint64_t)b[i] << (8*i);
    return v;
}


int main(void) {
    /* init constants */
    extern const uint64_t CONST_TABLE[];
    for (int i = 0; i < N_CONST; i++)
        value_vec[i] = CONST_TABLE[i];

    /* message = 8 bytes */
    uint8_t msg[8] = { 178,96,184,161,3,67,191,90 };

    /* expected digest */
    uint8_t digest_bytes[32] = {
        0xf3,0xa4,0x33,0x17,0x36,0xe3,0x9b,0x61,
        0xb5,0xea,0x9c,0x4f,0x07,0x66,0x12,0xe3,
        0xba,0xb7,0xd4,0xcf,0x36,0xce,0x5a,0xf1,
        0x86,0x31,0xd7,0xad,0xc5,0xde,0x69,0xe3
    };

    /* witness wires */
    Wire message = { N_CONST };
    Wire digest[4] = {
        { N_CONST+1 }, { N_CONST+2 },
        { N_CONST+3 }, { N_CONST+4 }
    };

    /* fill private inputs */
    value_vec[message.idx] = le_bytes_to_u64(msg, 8);
    for (int i = 0; i < 4; i++)
        value_vec[digest[i].idx] =
            le_bytes_to_u64(digest_bytes + 8*i, 8);

    /* padding (LEN=8 case) */
    Wire zero = { 0 };
    Wire padded[17];
    padded[0] = message;
    padded[1] = (Wire){1}; /* constant 0x01 */
    for (int i = 2; i < 17; i++) padded[i] = zero;
    padded[16] = (Wire){12}; /* 0x80<<56 */

    /* state init */
    Wire state[25];
    for (int i = 0; i < 25; i++) state[i] = zero;

    /* absorb */
    for (int i = 0; i < 17; i++)
        state[i] = (Wire){ next_internal++ }, value_vec[state[i].idx] =
            value_vec[padded[i].idx];

    /* permutation */
    keccak_f1600(state);

    /* check digest */
    for (int i = 0; i < 4; i++) {
        uint64_t got = value_vec[state[i].idx];
        uint64_t exp = value_vec[digest[i].idx];
        assert(got == exp);
    }

    printf("Digest OK ✔\n");

    /* dump witness */
    FILE *f = fopen("keccak_witness_dump_c.txt", "w");
    for (int i = 0; i < COMMITTED; i++)
        fprintf(f, "W[%05d] = 0x%016llx\n",
                i, (unsigned long long)value_vec[i]);
    fclose(f);

    printf("Witness written to keccak_witness_dump_c.txt\n");
    return 0;
}
