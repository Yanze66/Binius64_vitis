#include <ap_int.h>
#include <stdint.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include <ap_int.h>



using u8   = ap_uint<8>;
using u16  = ap_uint<16>;
using u32 = ap_uint<32>;
using u64  = ap_uint<64>;
using u128 = ap_uint<128>;
using u256 = ap_uint<256>;


typedef ap_axiu<32, 0, 0, 0>  axis32_t;
typedef ap_axiu<64, 0, 0, 0>  axis64_t;
typedef ap_axiu<128, 0, 0, 0> axis128_t;



static const int LOG_BITS = 6;
static const int HEIGHT   = 1 << LOG_BITS;   // 64
// static const int N_VARS   = 4;
static const int N_VARS   = 15;

static const int N        = 1 << N_VARS;     // 32768
static const int B_LEAVES_LEN = HEIGHT * N;  // 2097152

// ============================================================
// GF(2^128) GHASH arithmetic
// modulus: x^128 + x^7 + x^2 + x + 1
// ============================================================

static u64 reverse_bits_64(u64 x) {
#pragma HLS INLINE off//
    u64 r = 0;
    for (int i = 0; i < 64; i++) {
#pragma HLS UNROLL //不要用流水，就一个周期
// #pragma HLS PIPELINE II=1 // 流水，64个周期，面积更大，追求吞吐率
        r[63 - i] = x[i];
    }
    return r;
}


static u128 reverse_bits_each_64(u128 x) {
#pragma HLS INLINE off
    u64 lo = (u64)x;
    u64 hi = (u64)(x >> 64);

    //  #pragma HLS DATAFLOW //两个一起执行

    u64 lo_r = reverse_bits_64(lo);
    u64 hi_r = reverse_bits_64(hi);

    return (u128(hi_r) << 64) | u128(lo_r);
}

static u128 shr_each_64(u128 x, unsigned s) {
#pragma HLS INLINE off
    u64 lo = (u64)x;
    u64 hi = (u64)(x >> 64);

    
    lo >>= s;
    hi >>= s;

    return (u128(hi) << 64) | u128(lo);
}



// multiplication in GF(2), so no carry chain
// 流水线设计，每个clk迭代一次
static u128 clmul64(u64 a, u64 b) {
// #pragma HLS PIPELINE II=1

#pragma HLS INLINE off
    u128 acc = 0;
    for (int i = 0; i < 64; i++) {
        #pragma HLS UNROLL   // 完全展开，组合逻辑
        if (b[i]) {
            acc ^= (u128(a) << i);
        }
    }
    return acc;
}

// //块写法，尝试增加频率
// static u128 clmul64(u64 a, u64 b) {
// #pragma HLS INLINE off
//     u128 acc = 0;
//     for (int blk = 0; blk < 8; blk++) {
// #pragma HLS PIPELINE II=1
//         u128 part = 0;
//         for (int j = 0; j < 8; j++) {
// #pragma HLS UNROLL
//             int i = blk * 8 + j;
//             if (b[i]) {
//                 part ^= (u128(a) << i);
//             }
//         }
//         acc ^= part;
//     }
//     return acc;
// }


// 组合逻辑，看后续要不要变成pipeline，很有可能成为bottleneck
static u128 reduce_ghash_256_by_64(u64 v0, u64 v1, u64 v2, u64 v3) {
#pragma HLS INLINE 
// #pragma HLS PIPELINE II=1 //尝试流水线
    v1 ^= v3 ^ (v3 << 1) ^ (v3 << 2) ^ (v3 << 7);
    v2 ^= (v3 >> 63) ^ (v3 >> 62) ^ (v3 >> 57);
    v0 ^= v2 ^ (v2 << 1) ^ (v2 << 2) ^ (v2 << 7);
    v1 ^= (v2 >> 63) ^ (v2 >> 62) ^ (v2 >> 57);
    return (u128(v1) << 64) | u128(v0);
}

//GHASH irreducible polynomial: x^128 + x^7 + x^2 + x + 1 


/////////////////////////////////////////////////////////////
// pipelined mul ， improve its throughput 
/////////////////////////////////////////////////////////


static u128 ghash_mul_pipe_half(u128 x, u128 y) {
#pragma HLS INLINE 

    u64 x1 = (u64)(x >> 64);
    u64 x0 = (u64)(x);
    u64 y1 = (u64)(y >> 64);
    u64 y0 = (u64)(y);

    u64 x2 = x0 ^ x1;
    u64 y2 = y0 ^ y1;

    // only 3 carry-less multiplications
    u128 z0 = clmul64(y0, x0);
    u128 z1 = clmul64(y1, x1);
    u128 z2 = clmul64(y2, x2);

    // Karatsuba cross term
    z2 ^= z0 ^ z1;

    // assemble 256-bit product as 4x64-bit words
    u64 v0 = (u64)(z0);
    u64 v1 = (u64)(z0 >> 64) ^ (u64)(z2);
    u64 v2 = (u64)(z1) ^ (u64)(z2 >> 64);
    u64 v3 = (u64)(z1 >> 64);

    return reduce_ghash_256_by_64(v0, v1, v2, v3);
}



static u64 clmul32_base(u32 a, u32 b) {
#pragma HLS INLINE off
    u64 acc = 0;

    for (int blk = 0; blk < 8; blk++) {
#pragma HLS PIPELINE II=1
        u64 part = 0;

        for (int j = 0; j < 4; j++) {
#pragma HLS UNROLL
            int i = blk * 4 + j;
            if (b[i]) {
                part ^= (u64(a) << i);
            }
        }

        acc ^= part;
    }

    return acc;
}

static u128 clmul64_karatsuba32(u64 a, u64 b) {
#pragma HLS INLINE off

    u32 a0 = (u32)a;
    u32 a1 = (u32)(a >> 32);
    u32 b0 = (u32)b;
    u32 b1 = (u32)(b >> 32);

    u32 a2 = a0 ^ a1;
    u32 b2 = b0 ^ b1;

    u64 z0 = clmul32_base(a0, b0);
    u64 z1 = clmul32_base(a1, b1);
    u64 z2 = clmul32_base(a2, b2);

    z2 ^= z0 ^ z1;

    u128 res = 0;
    res ^= (u128)z0;
    res ^= ((u128)z2 << 32);
    res ^= ((u128)z1 << 64);

    return res;
}




static u128 clmul64_serial(u64 a, u64 b) {
#pragma HLS INLINE off
    u128 acc = 0;

    for (int blk = 0; blk < 16; blk++) {
#pragma HLS PIPELINE II=1
        u128 part = 0;

        for (int j = 0; j < 4; j++) {
#pragma HLS UNROLL
            int i = blk * 4 + j;
            if (b[i]) {
                part ^= (u128(a) << i);
            }
        }

        acc ^= part;
    }

    return acc;
}

static u128 ghash_mul_pipe_serial(u128 x, u128 y) {
#pragma HLS INLINE off

    u64 x1 = (u64)(x >> 64);
    u64 x0 = (u64)(x);
    u64 y1 = (u64)(y >> 64);
    u64 y0 = (u64)(y);

    u64 x2 = x0 ^ x1;
    u64 y2 = y0 ^ y1;

    // only 3 serial carry-less multiplications
    u128 z0 = clmul64_serial(y0, x0);
    u128 z1 = clmul64_serial(y1, x1);
    u128 z2 = clmul64_serial(y2, x2);
     

    // Karatsuba cross term
    z2 ^= z0 ^ z1;

    // assemble 256-bit product
    u64 v0 = (u64)(z0);
    u64 v1 = (u64)(z0 >> 64) ^ (u64)(z2);
    u64 v2 = (u64)(z1) ^ (u64)(z2 >> 64);
    u64 v3 = (u64)(z1 >> 64);

    return reduce_ghash_256_by_64(v0, v1, v2, v3);
}




static u128 ghash_mul_pipe_kara2(u128 x, u128 y) {
#pragma HLS INLINE off

    u64 x1 = (u64)(x >> 64);
    u64 x0 = (u64)(x);
    u64 y1 = (u64)(y >> 64);
    u64 y0 = (u64)(y);

    u64 x2 = x0 ^ x1;
    u64 y2 = y0 ^ y1;

    // only 3 serial carry-less multiplications
    
    u128 z0 = clmul64_karatsuba32(y0, x0);
    u128 z1 = clmul64_karatsuba32(y1, x1);
    u128 z2 = clmul64_karatsuba32(y2, x2);

    // Karatsuba cross term
    z2 ^= z0 ^ z1;

    // assemble 256-bit product
    u64 v0 = (u64)(z0);
    u64 v1 = (u64)(z0 >> 64) ^ (u64)(z2);
    u64 v2 = (u64)(z1) ^ (u64)(z2 >> 64);
    u64 v3 = (u64)(z1 >> 64);

    return reduce_ghash_256_by_64(v0, v1, v2, v3);
}

////////////////////////// 3-阶段karaatsuba ///////////////////////////////////////////////////


static u8 clmul4_lut(ap_uint<4> a, ap_uint<4> b) {
#pragma HLS INLINE off

    static const u8 table[256] = {
        // 预计算：clmul4(a, b)
        // index = (a << 4) | b

         // a = 0x0
        0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
        // a = 0x1
        0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
        // a = 0x2
        0x00,0x02,0x04,0x06,0x08,0x0a,0x0c,0x0e,0x10,0x12,0x14,0x16,0x18,0x1a,0x1c,0x1e,
        // a = 0x3
        0x00,0x03,0x06,0x05,0x0c,0x0f,0x0a,0x09,0x18,0x1b,0x1e,0x1d,0x14,0x17,0x12,0x11,

        // a = 0x4
        0x00,0x04,0x08,0x0c,0x10,0x14,0x18,0x1c,0x20,0x24,0x28,0x2c,0x30,0x34,0x38,0x3c,
        // a = 0x5
        0x00,0x05,0x0a,0x0f,0x14,0x11,0x1e,0x1b,0x28,0x2d,0x22,0x27,0x3c,0x39,0x36,0x33,
        // a = 0x6
        0x00,0x06,0x0c,0x0a,0x18,0x1e,0x14,0x12,0x30,0x36,0x3c,0x3a,0x28,0x2e,0x24,0x22,
        // a = 0x7
        0x00,0x07,0x0e,0x09,0x1c,0x1b,0x12,0x15,0x38,0x3f,0x36,0x31,0x24,0x23,0x2a,0x2d,

        // a = 0x8
        0x00,0x08,0x10,0x18,0x20,0x28,0x30,0x38,0x40,0x48,0x50,0x58,0x60,0x68,0x70,0x78,
        // a = 0x9
        0x00,0x09,0x12,0x1b,0x24,0x2d,0x36,0x3f,0x48,0x41,0x5a,0x53,0x6c,0x65,0x7e,0x77,
        // a = 0xa
        0x00,0x0a,0x14,0x1e,0x28,0x22,0x3c,0x36,0x50,0x5a,0x44,0x4e,0x78,0x72,0x6c,0x66,
        // a = 0xb
        0x00,0x0b,0x16,0x1d,0x2c,0x27,0x3a,0x31,0x58,0x53,0x4e,0x45,0x74,0x7f,0x62,0x69,

        // a = 0xc
        0x00,0x0c,0x18,0x14,0x30,0x3c,0x28,0x24,0x60,0x6c,0x78,0x74,0x50,0x5c,0x48,0x44,
        // a = 0xd
        0x00,0x0d,0x1a,0x17,0x34,0x39,0x2e,0x23,0x68,0x65,0x72,0x7f,0x5c,0x51,0x46,0x4b,
        // a = 0xe
        0x00,0x0e,0x1c,0x12,0x38,0x36,0x24,0x2a,0x70,0x7e,0x6c,0x62,0x48,0x46,0x54,0x5a,
        // a = 0xf
        0x00,0x0f,0x1e,0x11,0x3c,0x33,0x22,0x2d,0x78,0x77,0x66,0x69,0x44,0x4b,0x5a,0x55
    };

#pragma HLS BIND_STORAGE variable=table type=rom_1p impl=lutram

    ap_uint<8> idx = (a, b);
    return table[idx];
}

static u8 clmul4_comb(ap_uint<4> a, ap_uint<4> b) {
#pragma HLS INLINE off
    u8 acc = 0;

    for (int i = 0; i < 4; i++) {
#pragma HLS UNROLL
        if (b[i]) {
            acc ^= (u8(a) << i);
        }
    }

    return acc;
}


static u16 clmul8_partial4(u8 a, u8 b) {
#pragma HLS INLINE off
    u16 acc = 0;
    for (int blk = 0; blk < 2; blk++) {
#pragma HLS PIPELINE II=1
        u16 part = 0;
        for (int j = 0; j < 4; j++) {
#pragma HLS UNROLL
            int i = blk * 4 + j;
            if (b[i]) {
                part ^= (u16(a) << i);
            }
        }
        acc ^= part;
    }
    return acc;
}

static u16 clmul8_karatsuba4(u8 a, u8 b) {
#pragma HLS INLINE off

    ap_uint<4> a0 = a.range(3, 0);
    ap_uint<4> a1 = a.range(7, 4);
    ap_uint<4> b0 = b.range(3, 0);
    ap_uint<4> b1 = b.range(7, 4);

    ap_uint<4> a2 = a0 ^ a1;
    ap_uint<4> b2 = b0 ^ b1;

    u8 z0 = clmul4_lut(a0, b0);
    u8 z1 = clmul4_lut(a1, b1);
    u8 z2 = clmul4_lut(a2, b2);

    // u8 z0 = clmul4_comb(a0, b0);
    // u8 z1 = clmul4_comb(a1, b1);
    // u8 z2 = clmul4_comb(a2, b2);

    z2 ^= z0 ^ z1;

    u16 res = 0;
    res ^= (u16)z0;
    res ^= ((u16)z2 << 4);
    res ^= ((u16)z1 << 8);

    return res;
}

static u32 clmul16_karatsuba8(u16 a, u16 b) {
#pragma HLS INLINE off

    u8 a0 = (u8)a;
    u8 a1 = (u8)(a >> 8);
    u8 b0 = (u8)b;
    u8 b1 = (u8)(b >> 8);

    u8 a2 = a0 ^ a1;
    u8 b2 = b0 ^ b1;

    // u16 z0 = clmul8_partial4(a0, b0); 3阶段karatsuba
    // u16 z1 = clmul8_partial4(a1, b1);
    // u16 z2 = clmul8_partial4(a2, b2);
    u16 z0 = clmul8_karatsuba4(a0, b0); //4阶段karatsuba
    u16 z1 = clmul8_karatsuba4(a1, b1);
    u16 z2 = clmul8_karatsuba4(a2, b2);

    z2 ^= z0 ^ z1;

    u32 res = 0;
    res ^= (u32)z0;
    res ^= ((u32)z2 << 8);
    res ^= ((u32)z1 << 16);

    return res;
}

static u64 clmul32_karatsuba16(u32 a, u32 b) {
#pragma HLS INLINE off

    u16 a0 = (u16)a;
    u16 a1 = (u16)(a >> 16);
    u16 b0 = (u16)b;
    u16 b1 = (u16)(b >> 16);

    u16 a2 = a0 ^ a1;
    u16 b2 = b0 ^ b1;

    u32 z0 = clmul16_karatsuba8(a0, b0);
    u32 z1 = clmul16_karatsuba8(a1, b1);
    u32 z2 = clmul16_karatsuba8(a2, b2);

    z2 ^= z0 ^ z1;

    u64 res = 0;
    res ^= (u64)z0;
    res ^= ((u64)z2 << 16);
    res ^= ((u64)z1 << 32);

    return res;
}

static u128 clmul64_karatsuba(u64 a, u64 b) {
#pragma HLS INLINE off

    u32 a0 = (u32)a;
    u32 a1 = (u32)(a >> 32);
    u32 b0 = (u32)b;
    u32 b1 = (u32)(b >> 32);

    u32 a2 = a0 ^ a1;
    u32 b2 = b0 ^ b1;

    u64 z0 = clmul32_karatsuba16(a0, b0);
    u64 z1 = clmul32_karatsuba16(a1, b1);
    u64 z2 = clmul32_karatsuba16(a2, b2);

    z2 ^= z0 ^ z1;

    u128 res = 0;
    res ^= (u128)z0;
    res ^= ((u128)z2 << 32);
    res ^= ((u128)z1 << 64);

    return res;
}


static u128 ghash_mul_pipe_kara3(u128 x, u128 y) {
#pragma HLS INLINE off

    u64 x1 = (u64)(x >> 64);
    u64 x0 = (u64)(x);
    u64 y1 = (u64)(y >> 64);
    u64 y0 = (u64)(y);

    u64 x2 = x0 ^ x1;
    u64 y2 = y0 ^ y1;

    // only 3 serial carry-less multiplications
    
    u128 z0 = clmul64_karatsuba(y0, x0);
    u128 z1 = clmul64_karatsuba(y1, x1);
    u128 z2 = clmul64_karatsuba(y2, x2);

    // Karatsuba cross term
    z2 ^= z0 ^ z1;

    // assemble 256-bit product
    u64 v0 = (u64)(z0);
    u64 v1 = (u64)(z0 >> 64) ^ (u64)(z2);
    u64 v2 = (u64)(z1) ^ (u64)(z2 >> 64);
    u64 v3 = (u64)(z1 >> 64);

    return reduce_ghash_256_by_64(v0, v1, v2, v3);
}

////////////////////////// 4-阶段karaatsuba ///////////////////////////////////////////////////

static u128 ghash_mul_pipe_kara4(u128 x, u128 y) {
#pragma HLS INLINE off

    u64 x1 = (u64)(x >> 64);
    u64 x0 = (u64)(x);
    u64 y1 = (u64)(y >> 64);
    u64 y0 = (u64)(y);

    u64 x2 = x0 ^ x1;
    u64 y2 = y0 ^ y1;

    // only 3 serial carry-less multiplications
    
    u128 z0 = clmul64_karatsuba(y0, x0);
    u128 z1 = clmul64_karatsuba(y1, x1);
    u128 z2 = clmul64_karatsuba(y2, x2);

    // Karatsuba cross term
    z2 ^= z0 ^ z1;

    // assemble 256-bit product
    u64 v0 = (u64)(z0);
    u64 v1 = (u64)(z0 >> 64) ^ (u64)(z2);
    u64 v2 = (u64)(z1) ^ (u64)(z2 >> 64);
    u64 v3 = (u64)(z1 >> 64);

    return reduce_ghash_256_by_64(v0, v1, v2, v3);
}


////////////////////

static void write_axis128(hls::stream<axis128_t>& s, u128 x) {
#pragma HLS INLINE
    axis128_t v;
    v.data = (ap_uint<128>)x;
    s.write(v);
}



static u128 read_axis128(hls::stream<axis128_t>& s) {
#pragma HLS INLINE
    axis128_t v = s.read();
    return (u128)v.data;
}


// ============================================================
// Top function
// ============================================================
extern "C" {
void gf_mul_benchmark_top(
    hls::stream<axis128_t> &a_in,
    hls::stream<axis128_t> &b_in,
    hls::stream<axis128_t> &out
) {
#pragma HLS INTERFACE ap_ctrl_hs port=return
#pragma HLS INTERFACE axis port=a_in
#pragma HLS INTERFACE axis port=b_in
#pragma HLS INTERFACE axis port=out
#pragma HLS INTERFACE s_axilite port=return bundle=control

    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=4

        u128 a = read_axis128(a_in);
        u128 b = read_axis128(b_in);

        // u128 c = clmul256(a,b); // 原始组合逻辑 32798 cc, 40.1k lut , 136.2MHz
        // u128 c = ghash_mul_pipe(a, b); //32790 cyc， 46.8K LUT， 5.95 ns around 166MHz， 流水线优化，确保1个cyc处理
        // u128 c = ghash_mul_pipe_half(a, b); //    32790 cyc, 36.1K LUT, 194MHz  // 2026.4.1 综合面积比第一次小了很多，可能是后续kara实现中顺便优化掉了

        // u128 c = ghash_mul_pipe_serial_pipe(a, b); //98317 cyc ，16.7K LUT，189MHz, 去掉并路，面积换时间
        // u128 c = ghash_mul_pipe_kara2(a, b); //    32790 cyc, 28.1K LUT, 184MHz   2阶段kara，面积减少
        // u128 c = ghash_mul_pipe_kara3(a, b); // 32790 cyc, 20K LUT ， 158MHz， 3阶段kara，面积进一步减小
        u128 c = ghash_mul_pipe_kara4(a, b); // 32790 cyc, 18.6K LUT, 154MHz  4阶段kara     // karatsuba+查找表  32792,288MHz 20.1K LUT // 32834 cyc, 20.1K LUT , 579MHz //864MHz,12K FF, 9K LUT,  


        write_axis128(out, c);


    }
}
}