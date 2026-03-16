#include <ap_int.h>
#include <stdint.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include <ap_int.h>

using u64  = ap_uint<64>;
using u128 = ap_uint<128>;

typedef ap_axiu<64, 0, 0, 0>  axis64_t;
typedef ap_axiu<128, 0, 0, 0> axis128_t;

static const int LOG_BITS = 6;
static const int HEIGHT   = 1 << LOG_BITS;   // 64
static const int N_VARS   = 4;
// static const int N_VARS   = 15;

static const int N        = 1 << N_VARS;     // 32768
static const int B_LEAVES_LEN = HEIGHT * N;  // 2097152

// ============================================================
// GF(2^128) GHASH arithmetic
// modulus: x^128 + x^7 + x^2 + x + 1
// ============================================================

static u64 reverse_bits_64(u64 x) {
#pragma HLS INLINE
    u64 r = 0;
    for (int i = 0; i < 64; i++) {
#pragma HLS UNROLL
        r[63 - i] = x[i];
    }
    return r;
}


static u128 reverse_bits_each_64(u128 x) {
#pragma HLS INLINE
    u64 lo = (u64)x;
    u64 hi = (u64)(x >> 64);

    u64 lo_r = reverse_bits_64(lo);
    u64 hi_r = reverse_bits_64(hi);

    return (u128(hi_r) << 64) | u128(lo_r);
}

static u128 shr_each_64(u128 x, unsigned s) {
#pragma HLS INLINE
    u64 lo = (u64)x;
    u64 hi = (u64)(x >> 64);

    lo >>= s;
    hi >>= s;

    return (u128(hi) << 64) | u128(lo);
}

//multiplication in GF(2), so no carry chain
static u128 clmul64(u64 a, u64 b) {
#pragma HLS INLINE
    u128 acc = 0;
    for (int i = 0; i < 64; i++) {
#pragma HLS UNROLL factor=1
        if (b[i]) {
            acc ^= (u128(a) << i);
        }
    }
    return acc;
}

static u128 reduce_ghash_256_by_64(u64 v0, u64 v1, u64 v2, u64 v3) {
#pragma HLS INLINE
    v1 ^= v3 ^ (v3 << 1) ^ (v3 << 2) ^ (v3 << 7);
    v2 ^= (v3 >> 63) ^ (v3 >> 62) ^ (v3 >> 57);
    v0 ^= v2 ^ (v2 << 1) ^ (v2 << 2) ^ (v2 << 7);
    v1 ^= (v2 >> 63) ^ (v2 >> 62) ^ (v2 >> 57);
    return (u128(v1) << 64) | u128(v0);
}

static u128 ghash_mul(u128 x, u128 y) {
#pragma HLS INLINE
    u64 x1 = (u64)(x >> 64);
    u64 x0 = (u64)(x);
    u64 y1 = (u64)(y >> 64);
    u64 y0 = (u64)(y);

    u64 x0r = reverse_bits_64(x0);
    u64 x1r = reverse_bits_64(x1);
    u64 x2  = x0 ^ x1;

    u64 y0r = reverse_bits_64(y0);
    u64 y1r = reverse_bits_64(y1);
    u64 y2  = y0 ^ y1;

    u128 z0  = clmul64(y0,  x0);
    u128 z1  = clmul64(y1,  x1);
    u128 z2  = clmul64(y2,  x2);

    u128 z0h = clmul64(y0r, x0r);
    u128 z1h = clmul64(y1r, x1r);
    u128 z2h = clmul64(y0r ^ y1r, x0r ^ x1r);

    z2  ^= z0 ^ z1;
    z2h ^= z0h ^ z1h;

z0h = shr_each_64(reverse_bits_each_64(z0h), 1);
z1h = shr_each_64(reverse_bits_each_64(z1h), 1);
z2h = shr_each_64(reverse_bits_each_64(z2h), 1);

    u64 v0 = (u64)(z0);
    u64 v1 = (u64)(z0h) ^ (u64)(z2);
    u64 v2 = (u64)(z1)  ^ (u64)(z2h);
    u64 v3 = (u64)(z1h);

    return reduce_ghash_256_by_64(v0, v1, v2, v3);
}

static inline u128 gf_add(u128 a, u128 b) {
#pragma HLS INLINE
    return a ^ b;
}

static inline u128 gf_square(u128 a) {
#pragma HLS INLINE
    return ghash_mul(a, a);
}

static inline u128 gf_one() {
#pragma HLS INLINE
    return (u128)1;
}

static u128 gf_pow_u64(u128 base, u64 exp) {
#pragma HLS INLINE off
    u128 result = gf_one();
    u128 cur = base;
    for (int i = 0; i < 64; i++) {
#pragma HLS PIPELINE II=1
        if (exp[i]) {
            result = ghash_mul(result, cur);
        }
        cur = gf_square(cur);
    }
    return result;
}

// ============================================================
// GHASH field constants from Rust
// ============================================================

// F::MULTIPLICATIVE_GENERATOR for BinaryField128bGhash
static const u128 GHASH_GENERATOR =
    (((u128)0x494ef99794d5244full) << 64) |
     (u128)0x9152df59d87a9186ull;

// g_c_hi = iterate(g, |g| g.square()).nth(1 << log_bits)
// log_bits = 6 => 64 squarings
static u128 compute_g_c_hi() {
#pragma HLS INLINE off
    u128 g = GHASH_GENERATOR;
    for (int i = 0; i < HEIGHT; i++) {
#pragma HLS PIPELINE II=1
        g = gf_square(g);
    }
    return g;
}

// ============================================================
// Build root table for BinaryTree::constant_base(log_bits=6, base, exponents)
// Since the full tree root equals base^(exponent[i]) per vertex, we can compute
// the final root directly.
// ============================================================

static void build_constant_base_root(
    const u64 exponents[N],
    u128 base,
    u128 root_out[N]
) {
#pragma HLS INLINE off
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        root_out[i] = gf_pow_u64(base, exponents[i]);
    }
}

// ============================================================
// compute_b_leaves(log_bits, bases=a_root, exponents=b)
// out[z * N + i] = bit z of b[i] ? (a_root[i])^(2^z) : 1
// ============================================================

static void build_b_leaves(
    const u128 a_root[N],
    const u64 b_exp[N],
    u128 b_leaves[B_LEAVES_LEN]
) {
#pragma HLS INLINE off
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE off
        u128 cur = a_root[i];
        u64 exp  = b_exp[i];

        for (int z = 0; z < HEIGHT; z++) {
#pragma HLS PIPELINE II=1
            bool bit = (bool)((exp >> z) & 1);
            b_leaves[z * N + i] = bit ? cur : gf_one();
            cur = gf_square(cur);
        }
    }
}

// ============================================================
// ProdcheckProver::new(k=6, witness=b_leaves)
// Final products layer => b_root
//
// witness length = 2^(n_vars + log_bits) = 2^21
// each reduction halves the front dimension by pairwise product of split halves
// after 6 rounds => length 2^15 = 32768
// ============================================================

// static void build_b_root_from_b_leaves(
//     const u128 b_leaves[B_LEAVES_LEN],
//     u128 b_root[N]
// ) {
// #pragma HLS INLINE off

//     static u128 layer0[B_LEAVES_LEN];
//     static u128 layer1[B_LEAVES_LEN / 2];
//     static u128 layer2[B_LEAVES_LEN / 4];
//     static u128 layer3[B_LEAVES_LEN / 8];
//     static u128 layer4[B_LEAVES_LEN / 16];
//     static u128 layer5[B_LEAVES_LEN / 32];
// #pragma HLS BIND_STORAGE variable=layer0 type=ram_2p impl=uram
// #pragma HLS BIND_STORAGE variable=layer1 type=ram_2p impl=uram
// #pragma HLS BIND_STORAGE variable=layer2 type=ram_2p impl=uram
// #pragma HLS BIND_STORAGE variable=layer3 type=ram_2p impl=uram
// #pragma HLS BIND_STORAGE variable=layer4 type=ram_2p impl=uram
// #pragma HLS BIND_STORAGE variable=layer5 type=ram_2p impl=uram

//     for (int i = 0; i < B_LEAVES_LEN; i++) {
// #pragma HLS PIPELINE II=1
//         layer0[i] = b_leaves[i];
//     }

//     // 2^21 -> 2^20
//     for (int i = 0; i < (B_LEAVES_LEN / 2); i++) {
// #pragma HLS PIPELINE II=1
//         layer1[i] = ghash_mul(layer0[i], layer0[i + (B_LEAVES_LEN / 2)]);
//     }

//     // 2^20 -> 2^19
//     for (int i = 0; i < (B_LEAVES_LEN / 4); i++) {
// #pragma HLS PIPELINE II=1
//         layer2[i] = ghash_mul(layer1[i], layer1[i + (B_LEAVES_LEN / 4)]);
//     }

//     // 2^19 -> 2^18
//     for (int i = 0; i < (B_LEAVES_LEN / 8); i++) {
// #pragma HLS PIPELINE II=1
//         layer3[i] = ghash_mul(layer2[i], layer2[i + (B_LEAVES_LEN / 8)]);
//     }

//     // 2^18 -> 2^17
//     for (int i = 0; i < (B_LEAVES_LEN / 16); i++) {
// #pragma HLS PIPELINE II=1
//         layer4[i] = ghash_mul(layer3[i], layer3[i + (B_LEAVES_LEN / 16)]);
//     }

//     // 2^17 -> 2^16
//     for (int i = 0; i < (B_LEAVES_LEN / 32); i++) {
// #pragma HLS PIPELINE II=1
//         layer5[i] = ghash_mul(layer4[i], layer4[i + (B_LEAVES_LEN / 32)]);
//     }

//     // 2^16 -> 2^15 = b_root
//     for (int i = 0; i < N; i++) {
// #pragma HLS PIPELINE II=1
//         b_root[i] = ghash_mul(layer5[i], layer5[i + N]);
//     }
// }

// ============================================================
// c_root = c_lo.root * c_hi.root
// ============================================================

static void build_c_root(
    const u128 clo_root[N],
    const u128 chi_root[N],
    u128 c_root[N]
) {
#pragma HLS INLINE off
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        c_root[i] = ghash_mul(clo_root[i], chi_root[i]);
    }
}


static u64 read_axis64(hls::stream<axis64_t>& s) {
#pragma HLS INLINE
    axis64_t v = s.read();
    return (u64)v.data;
}

static void write_axis128(hls::stream<axis128_t>& s, u128 x) {
#pragma HLS INLINE
    axis128_t v;
    v.data = (ap_uint<128>)x;
    s.write(v);
}

// ============================================================
// Top function
// ============================================================

// extern "C" {
// void intmul_witness_step(
//     hls::stream<axis64_t> &a_in,
//     hls::stream<axis64_t> &b_in,
//     hls::stream<axis64_t> &clo_in,
//     hls::stream<axis64_t> &chi_in,

//     hls::stream<axis128_t> &a_root_out,
//     hls::stream<axis128_t> &clo_root_out,
//     hls::stream<axis128_t> &chi_root_out,
//     hls::stream<axis128_t> &b_leaves_out
//     // hls::stream<axis128_t> &b_root_out,
//     // hls::stream<axis128_t> &c_root_out
// ) {
// #pragma HLS INTERFACE ap_ctrl_hs port=return

// #pragma HLS INTERFACE axis port=a_in
// #pragma HLS INTERFACE axis port=b_in
// #pragma HLS INTERFACE axis port=clo_in
// #pragma HLS INTERFACE axis port=chi_in

// #pragma HLS INTERFACE axis port=a_root_out
// #pragma HLS INTERFACE axis port=clo_root_out
// #pragma HLS INTERFACE axis port=chi_root_out
// #pragma HLS INTERFACE axis port=b_leaves_out
// // #pragma HLS INTERFACE axis port=b_root_out
// // #pragma HLS INTERFACE axis port=c_root_out

// #pragma HLS INTERFACE s_axilite port=return bundle=control

//     static u64 a_raw[N];
//     static u64 b_raw[N];
//     static u64 clo_raw[N];
//     static u64 chi_raw[N];

//     static u128 a_root[N];
//     static u128 clo_root[N];
//     static u128 chi_root[N];
//     static u128 b_leaves[B_LEAVES_LEN];
//     // static u128 b_root[N];
//     // static u128 c_root[N];

// #pragma HLS BIND_STORAGE variable=a_raw type=ram_2p impl=bram
// #pragma HLS BIND_STORAGE variable=b_raw type=ram_2p impl=bram
// #pragma HLS BIND_STORAGE variable=clo_raw type=ram_2p impl=bram
// #pragma HLS BIND_STORAGE variable=chi_raw type=ram_2p impl=bram

// #pragma HLS BIND_STORAGE variable=a_root type=ram_2p impl=uram
// #pragma HLS BIND_STORAGE variable=clo_root type=ram_2p impl=uram
// #pragma HLS BIND_STORAGE variable=chi_root type=ram_2p impl=uram
// #pragma HLS BIND_STORAGE variable=b_leaves type=ram_2p impl=uram
// // #pragma HLS BIND_STORAGE variable=b_root type=ram_2p impl=uram
// // #pragma HLS BIND_STORAGE variable=c_root type=ram_2p impl=uram

//     // -------------------------
//     // Read 4 raw witness tables
//     // -------------------------
//     for (int i = 0; i < N; i++) {
// #pragma HLS PIPELINE II=1
//         a_raw[i]   = read_axis64(a_in);
//         b_raw[i]   = read_axis64(b_in);
//         clo_raw[i] = read_axis64(clo_in);
//         chi_raw[i] = read_axis64(chi_in);
//     }

//     u128 g = GHASH_GENERATOR;
//     u128 g_c_hi = compute_g_c_hi();

//     // a.root = g^(a[i])
//     build_constant_base_root(a_raw, g, a_root);

//     // c_lo.root = g^(c_lo[i])
//     build_constant_base_root(clo_raw, g, clo_root);

//     // c_hi.root = g_c_hi^(c_hi[i])
//     build_constant_base_root(chi_raw, g_c_hi, chi_root);

//     // b_leaves from a_root and b
//     build_b_leaves(a_root, b_raw, b_leaves);

//     // b_root by prodcheck reduction over b_leaves
//     // build_b_root_from_b_leaves(b_leaves, b_root);

//     // c_root = c_lo.root * c_hi.root
//     // build_c_root(clo_root, chi_root, c_root);

//     // -------------------------
//     // Stream out results
//     // -------------------------
//     for (int i = 0; i < N; i++) {
// #pragma HLS PIPELINE II=1
//         write_axis128(a_root_out, a_root[i]);
//         write_axis128(clo_root_out, clo_root[i]);
//         write_axis128(chi_root_out, chi_root[i]);
//     }

//     for (int i = 0; i < B_LEAVES_LEN; i++) {
// #pragma HLS PIPELINE II=1
//         write_axis128(b_leaves_out, b_leaves[i]);
//     }

// //     for (int i = 0; i < N; i++) {
// // #pragma HLS PIPELINE II=1
// //         write_axis128(b_root_out, b_root[i]);
// //         write_axis128(c_root_out, c_root[i]);
// //     }
// }
// }

// systhesis get stuck, so I simplize the kernel

extern "C" {
void intmul_witness_step(
    hls::stream<axis64_t> &a_in,
    hls::stream<axis64_t> &b_in,
    hls::stream<axis128_t> &a_root_out,
    hls::stream<axis128_t> &b_leaves_out
) {
#pragma HLS INTERFACE ap_ctrl_hs port=return
#pragma HLS INTERFACE axis port=a_in
#pragma HLS INTERFACE axis port=b_in
#pragma HLS INTERFACE axis port=a_root_out
#pragma HLS INTERFACE axis port=b_leaves_out
#pragma HLS INTERFACE s_axilite port=return bundle=control

    static u64 a_raw[N];
    static u64 b_raw[N];
    static u128 a_root[N];
    static u128 b_leaves[B_LEAVES_LEN];

#pragma HLS BIND_STORAGE variable=a_raw type=ram_2p impl=bram
#pragma HLS BIND_STORAGE variable=b_raw type=ram_2p impl=bram
#pragma HLS BIND_STORAGE variable=a_root type=ram_2p impl=bram
#pragma HLS BIND_STORAGE variable=b_leaves type=ram_2p impl=uram

    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        a_raw[i] = read_axis64(a_in);
        b_raw[i] = read_axis64(b_in);
    }

    build_constant_base_root(a_raw, GHASH_GENERATOR, a_root);
    build_b_leaves(a_root, b_raw, b_leaves);

    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        write_axis128(a_root_out, a_root[i]);
    }

    for (int i = 0; i < B_LEAVES_LEN; i++) {
#pragma HLS PIPELINE II=1
        write_axis128(b_leaves_out, b_leaves[i]);
    }
}
}
