#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* ---------- constants (same as Rust) ---------- */

static const uint64_t RC[24] = {
    0x0000000000000001ULL,
    0x0000000000008082ULL,
    0x800000000000808aULL,
    0x8000000080008000ULL,
    0x000000000000808bULL,
    0x0000000080000001ULL,
    0x8000000080008081ULL,
    0x8000000000008009ULL,
    0x000000000000008aULL,
    0x0000000000000088ULL,
    0x0000000080008009ULL,
    0x000000008000000aULL,
    0x000000008000808bULL,
    0x800000000000008bULL,
    0x8000000000008089ULL,
    0x8000000000008003ULL,
    0x8000000000008002ULL,
    0x8000000000000080ULL,
    0x000000000000800aULL,
    0x800000008000000aULL,
    0x8000000080008081ULL,
    0x8000000000008080ULL,
    0x0000000080000001ULL,
    0x8000000080008008ULL,
};

static const int R[25] = {
     0,  1, 62, 28, 27,
    36, 44,  6, 55, 20,
     3, 10, 43, 25, 39,
    41, 45, 15, 21,  8,
    18,  2, 61, 56, 14
};

static inline int idx(int x, int y) {
    return x + 5 * y;
}

static inline uint64_t rotl(uint64_t x, int n) {
    return (x << n) | (x >> (64 - n));
}

static void dump_a0(const char* tag, uint64_t A[25]) {
    printf("%s A[0]=0x%016llx\n", tag, (unsigned long long)A[0]);
}




/* ---------- keccak steps ---------- */

static void theta(uint64_t A[25]) {
    uint64_t C[5], D[5];

    for (int x = 0; x < 5; x++) {
        C[x] = A[idx(x,0)] ^ A[idx(x,1)] ^ A[idx(x,2)]
             ^ A[idx(x,3)] ^ A[idx(x,4)];
    }

    for (int x = 0; x < 5; x++) {
        D[x] = C[(x+4)%5] ^ rotl(C[(x+1)%5], 1);
    }

	/* ===== Rust force_commit 对应点 ===== */
    printf("=== ROUND %d : THETA force_commit D[x] ===\n");
    for (int x = 0; x < 5; x++) {
        printf("D[%d] = 0x%016llx\n",
               x, (unsigned long long)D[x]);
    }

    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {
            A[idx(x,y)] ^= D[x];
        }
    }
}

static void rho_pi(uint64_t A[25]) {
    uint64_t T[25];
    memcpy(T, A, sizeof(T));

    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {
            int i = idx(x,y);
            int j = idx(y, (2*x + 3*y) % 5);
            if (R[i] == 0)
                T[j] = A[i];
            else
                T[j] = rotl(A[i], R[i]);
        }
    }

    memcpy(A, T, sizeof(T));
}

static void chi(uint64_t A[25]) {
    uint64_t T[5];
    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++)
            T[x] = A[idx(x,y)];
        for (int x = 0; x < 5; x++)
            A[idx(x,y)] ^= (~T[(x+1)%5]) & T[(x+2)%5];
    }
}

static void iota(uint64_t A[25], int round) {
    A[0] ^= RC[round];
}

static void keccak_round(uint64_t A[25], int round) {
    //dump_a0("before theta:", A);
    theta(A);
    //dump_a0("before rho:", A);
    rho_pi(A);
    //dump_a0("before chi:", A);
    chi(A);
    //dump_a0("before theta:", A);
    /* dump state */
    printf("=== STATE AFTER chi ===\n");
    for (int i = 0; i < 25; i++) {
        printf("A[%02d] = 0x%016llx\n",
               i, (unsigned long long)A[i]);
    }
    iota(A, round);
    //dump_a0("final A[0]:", A);
}

size_t pack_message_le(
    const uint8_t *msg,
    size_t len_bytes,
    uint64_t *out_words
) {
    size_t n_words = (len_bytes + 7) / 8;
    memset(out_words, 0, n_words * sizeof(uint64_t));

    for (size_t i = 0; i < len_bytes; i++) {
        size_t w = i / 8;
        size_t b = i % 8;
        out_words[w] |= ((uint64_t)msg[i]) << (8 * b);
    }
    return n_words;
}

/* ---------- test: EXACTLY like Rust ---------- */

int main(void) {
    uint64_t state[25];
    memset(state, 0, sizeof(state));

    /* ================= input ================= */

    /* ===== choose message ===== */

    // --- 1 byte ---
    //uint8_t message[] = {0x61};
    //size_t len_bytes = 1;

    // --- 8 bytes ---
    uint8_t message[] = {0xb2,0x60,0xb8,0xa1,0x03,0x43,0xbf,0x5a};
    size_t len_bytes = 8;

    // --- 135 bytes  ---
    //uint8_t message[135] = {
    //    205,56,46,120,38,70,176,74, 76,161,62,248,186,171,0,97,
    //    14,97,19,181,18,218,255,110, 190,82,13,172,78,15,252,212,
    //    167,160,145,109,245,59,123,1, 104,252,210,46,223,161,102,93,
    //    148,109,214,57,223,206,145,171, 245,78,42,68,87,119,185,75,
    //    101,248,209,146,155,53,143,129, 173,242,54,174,209,171,27,215,
    //    198,250,92,195,67,161,79,113, 156,60,94,138,44,79,95,120,
    //    79,40,68,247,42,43,19,68, 34,230,149,237,238,160,33,80,
    //    23,97,248,116,231,55,174,141, 240,103,50,221,30,28,239,242,
    //    231,101,207,149,236,37,35
    //};
    //size_t len_bytes = 135;

    /* ================= padding ================= */
    const size_t RATE_WORDS = 17;

    /* number of blocks */
    size_t n_blocks = (len_bytes + 1 + 136 - 1) / 136;
    size_t n_padded_words = n_blocks * RATE_WORDS;

    /* padded message */
    uint64_t *padded = calloc(n_padded_words, sizeof(uint64_t));


        size_t full_words = len_bytes / 8;
    size_t rem_bytes  = len_bytes % 8;

    /* full words */
    for (size_t i = 0; i < full_words; i++) {
        uint64_t w = 0;
        for (int b = 0; b < 8; b++)
            w |= ((uint64_t)message[8*i + b]) << (8*b);
        padded[i] = w;
    }

    /* boundary word */
    if (rem_bytes == 0) {
        padded[full_words] = 0x01ULL;
    } else {
        uint64_t last = 0;
        for (size_t b = 0; b < rem_bytes; b++)
            last |= ((uint64_t)message[8*full_words + b]) << (8*b);

        uint64_t padding_bit = 1ULL << (rem_bytes * 8);
        padded[full_words] = last ^ padding_bit;
    }

    /* ---- final 0x80 << 56 ---- */
    padded[n_blocks * RATE_WORDS - 1] ^= (0x80ULL << 56);

    /* ================= debug padding ================= */
    //printf("=== PADDED WORDS ===\n");
    //for (int i = 0; i < RATE_WORDS; i++) {
    //    printf("padded[%02d] = 0x%016llx\n",
    //           i, (unsigned long long)padded[i]);
    //}

    /* ================= absorb and permutation ================= */
    for (size_t b = 0; b < n_blocks; b++) {
    /* absorb one block */
    for (int i = 0; i < RATE_WORDS; i++) {
        state[i] ^= padded[b * RATE_WORDS + i];
    }

    /* permutation */
    for (int r = 0; r < 24; r++) {
        keccak_round(state, r);
    }
    }

    /* ================= digest ================= */
    printf("\n=== DIGEST (state[0..3]) ===\n");
    for (int i = 0; i < 4; i++) {
        printf("digest[%d] = 0x%016llx\n",
               i, (unsigned long long)state[i]);
    }

    return 0;
}
