#include <cstdint>
#include <vector>
#include <cstdio>
#include <cstring>

static inline uint64_t rotl64(uint64_t x, uint32_t r) {
    return (x << r) | (x >> (64 - r));
}

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
    0x0000000080000001ULL, 0x8000000080008008ULL,
};

static const uint32_t R[25] = {
     0,  1, 62, 28, 27,
    36, 44,  6, 55, 20,
     3, 10, 43, 25, 39,
    41, 45, 15, 21,  8,
    18,  2, 61, 56, 14
};

#define IDX(x,y) ((x) + 5*(y))

int main() {
    std::vector<uint64_t> W;
    W.reserve(1024);

    /* === constants (25) === */
    for (int i = 0; i < 25; i++)
        W.push_back(0);

    /* === input witness === */
    uint64_t message = 0x5abf4303a1b860b2ULL;
    W.push_back(message);

    /* === expected digest witness (copied from Rust test output) === */
    W.push_back(0x619be3361733a4f3ULL);
    W.push_back(0xe31266074f9ceab5ULL);
    W.push_back(0xf15ace36cfd4b7baULL);
    W.push_back(0xe369dec5add73186ULL);

    /* === internal trace === */
    uint64_t A[25] = {0};
    A[0] ^= message;

    for (int round = 0; round < 24; round++) {

        /* ---- theta ---- */
        uint64_t C[5], D[5];

        for (int x = 0; x < 5; x++) {
            C[x] = A[IDX(x,0)] ^ A[IDX(x,1)] ^ A[IDX(x,2)]
                 ^ A[IDX(x,3)] ^ A[IDX(x,4)];
            W.push_back(C[x]);  // bxor_multi
        }

        for (int x = 0; x < 5; x++) {
            uint64_t rot = rotl64(C[(x + 1) % 5], 1);
            W.push_back(rot);   // rotl
            D[x] = C[(x + 4) % 5] ^ rot;
            W.push_back(D[x]);  // bxor
        }

        for (int y = 0; y < 5; y++)
            for (int x = 0; x < 5; x++) {
                A[IDX(x,y)] ^= D[x];
                W.push_back(A[IDX(x,y)]);
            }

        /* ---- rho + pi ---- */
        uint64_t B[25];
        memcpy(B, A, sizeof(A));

        for (int y = 0; y < 5; y++)
            for (int x = 0; x < 5; x++) {
                uint32_t r = R[IDX(x,y)];
                if (r != 0) {
                    uint64_t v = rotl64(A[IDX(x,y)], r);
                    B[IDX(y, (2*x + 3*y) % 5)] = v;
                    W.push_back(v);
                }
            }
        memcpy(A, B, sizeof(A));

        /* ---- chi ---- */
        for (int y = 0; y < 5; y++) {
            uint64_t row[5];
            for (int x = 0; x < 5; x++)
                row[x] = A[IDX(x,y)];

            for (int x = 0; x < 5; x++) {
                uint64_t notb = ~row[(x + 1) % 5];
                W.push_back(notb);
                uint64_t andv = notb & row[(x + 2) % 5];
                W.push_back(andv);
                A[IDX(x,y)] = row[x] ^ andv;
                W.push_back(A[IDX(x,y)]);
            }
        }

        /* ---- iota ---- */
        A[0] ^= RC[round];
        W.push_back(A[0]);
    }

    /* === pad to committed size === */
    while (W.size() < 1024)
        W.push_back(0);

    /* === dump === */
    FILE* f = fopen("keccak_witness_dump.txt", "w");
    for (size_t i = 0; i < W.size(); i++) {
        fprintf(f, "W[%04zu] = 0x%016lx\n", i, W[i]);
    }
    fclose(f);

    return 0;
}