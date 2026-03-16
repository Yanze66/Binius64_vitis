#include <ap_int.h>
#include <ap_axi_sdata.h>
#include <hls_stream.h>

#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <limits.h>

using u64  = ap_uint<64>;
using u128 = ap_uint<128>;

typedef ap_axiu<64, 0, 0, 0>  axis64_t;
typedef ap_axiu<128, 0, 0, 0> axis128_t;

static const int LOG_BITS = 6;
static const int HEIGHT   = 1 << LOG_BITS;
static const int N_VARS   = 4;
// static const int N_VARS   = 15;


// 调试时你可以先改成 16
static const int N = 1 << N_VARS;

static const int B_LEAVES_LEN = HEIGHT * N;

extern "C" void intmul_witness_step(
    hls::stream<axis64_t> &a_in,
    hls::stream<axis64_t> &b_in,
    // hls::stream<axis64_t> &clo_in,
    // hls::stream<axis64_t> &chi_in,

    hls::stream<axis128_t> &a_root_out,
    // hls::stream<axis128_t> &clo_root_out,
    // hls::stream<axis128_t> &chi_root_out,
    hls::stream<axis128_t> &b_leaves_out
    // hls::stream<axis128_t> &b_root_out,
    // hls::stream<axis128_t> &c_root_out
);

static const char* A_FILE  = "intmul_witness_a.txt";
static const char* B_FILE  = "intmul_witness_b.txt";
static const char* LO_FILE = "intmul_witness_clo.txt";
static const char* HI_FILE = "intmul_witness_chi.txt";

static void print_cwd() {
    char buf[PATH_MAX];
    if (getcwd(buf, sizeof(buf)) != nullptr) {
        std::cout << "[TB] current working directory = " << buf << "\n";
    }
}

static std::string trim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace((unsigned char)s[l])) l++;
    size_t r = s.size();
    while (r > l && std::isspace((unsigned char)s[r - 1])) r--;
    return s.substr(l, r - l);
}

static bool starts_with_hex_index_line(const std::string& line) {
    if (line.empty()) return false;
    size_t i = 0;
    while (i < line.size() && std::isdigit((unsigned char)line[i])) i++;
    if (i == 0) return false;
    while (i < line.size() && std::isspace((unsigned char)line[i])) i++;
    return (i + 2 <= line.size() && line[i] == '0' && (line[i + 1] == 'x' || line[i + 1] == 'X'));
}

static bool parse_index_hex_u64_line(const std::string& line, int& idx, uint64_t& value) {
    std::stringstream ss(line);
    std::string idx_str, hex_str;
    if (!(ss >> idx_str >> hex_str)) return false;

    try { idx = std::stoi(idx_str); } catch (...) { return false; }

    if (!(hex_str.size() >= 3 && hex_str[0] == '0' && (hex_str[1] == 'x' || hex_str[1] == 'X')))
        return false;

    try { value = std::stoull(hex_str, nullptr, 16); } catch (...) { return false; }
    return true;
}

static bool load_witness_u64_file(const std::string& path, u64 out[N]) {
    std::ifstream fin(path);
    if (!fin) {
        std::cerr << "ERROR: cannot open file: " << path << "\n";
        return false;
    }

    std::string line;
    int count = 0;
    while (std::getline(fin, line)) {
        line = trim(line);
        if (!starts_with_hex_index_line(line)) continue;

        int idx_dummy = -1;
        uint64_t val = 0;
        if (!parse_index_hex_u64_line(line, idx_dummy, val)) continue;

        out[count] = (u64)val;
        count++;

        if (count == N) break;
    }

    if (count < N) {
        std::cerr << "ERROR: file " << path
                  << " only contains " << count
                  << " values, expected at least " << N << "\n";
        return false;
    }

    return true;
}

static std::string u64_to_hex(u64 x) {
    uint64_t lo = (uint64_t)x;
    std::ostringstream oss;
    oss << "0x" << std::hex << std::setfill('0') << std::setw(16) << lo;
    return oss.str();
}

static std::string u128_to_hex(u128 x) {
    uint64_t hi = (uint64_t)(x >> 64);
    uint64_t lo = (uint64_t)(x);
    std::ostringstream oss;
    oss << "0x"
        << std::hex << std::setfill('0')
        << std::setw(16) << hi
        << std::setw(16) << lo;
    return oss.str();
}

static void push_axis64(hls::stream<axis64_t>& s, u64 x) {
    axis64_t v;
    v.data = (ap_uint<64>)x;
    s.write(v);
}

static u128 pop_axis128(hls::stream<axis128_t>& s) {
    axis128_t v = s.read();
    return (u128)v.data;
}

int main() {
    print_cwd();

    static u64 a_raw[N];
    static u64 b_raw[N];
    // static u64 clo_raw[N];
    // static u64 chi_raw[N];

    static u128 a_root[N];
    // static u128 clo_root[N];
    // static u128 chi_root[N];
    static u128 b_leaves[B_LEAVES_LEN];
    // static u128 b_root[N];
    // static u128 c_root[N];

    if (!load_witness_u64_file(A_FILE, a_raw))   return 1;
    if (!load_witness_u64_file(B_FILE, b_raw))   return 1;
    // if (!load_witness_u64_file(LO_FILE, clo_raw)) return 1;
    // if (!load_witness_u64_file(HI_FILE, chi_raw)) return 1;

    hls::stream<axis64_t> a_in("a_in");
    hls::stream<axis64_t> b_in("b_in");
    // hls::stream<axis64_t> clo_in("clo_in");
    // hls::stream<axis64_t> chi_in("chi_in");

    hls::stream<axis128_t> a_root_out("a_root_out");
    // hls::stream<axis128_t> clo_root_out("clo_root_out");
    // hls::stream<axis128_t> chi_root_out("chi_root_out");
    hls::stream<axis128_t> b_leaves_out("b_leaves_out");
    // hls::stream<axis128_t> b_root_out("b_root_out");
    // hls::stream<axis128_t> c_root_out("c_root_out");

    for (int i = 0; i < N; i++) {
        push_axis64(a_in,   a_raw[i] );
        push_axis64(b_in,   b_raw[i] );
        // push_axis64(clo_in, clo_raw[i]);
        // push_axis64(chi_in, chi_raw[i]);
    }

    // intmul_witness_step(
    //     a_in, b_in, clo_in, chi_in,
    //     a_root_out, clo_root_out, chi_root_out,
    //     b_leaves_out, b_root_out, c_root_out
    // );
    intmul_witness_step(a_in, b_in, a_root_out, b_leaves_out);
    
    for (int i = 0; i < N; i++) a_root[i] = pop_axis128(a_root_out);
    // for (int i = 0; i < N; i++) clo_root[i] = pop_axis128(clo_root_out);
    // for (int i = 0; i < N; i++) chi_root[i] = pop_axis128(chi_root_out);
    for (int i = 0; i < B_LEAVES_LEN; i++) b_leaves[i] = pop_axis128(b_leaves_out);
    // for (int i = 0; i < N; i++) b_root[i] = pop_axis128(b_root_out);
    // for (int i = 0; i < N; i++) c_root[i] = pop_axis128(c_root_out);

    std::cout << "\n== a_root first 10 ==\n";
    for (int i = 0; i < (N < 10 ? N : 10); i++)
        std::cout << i << "\t" << u128_to_hex(a_root[i]) << "\n";

    std::cout << "\n== b-leaves first 10 ==\n";
    for (int i = 0; i < (N < 10 ? N : 10); i++)
        std::cout << i << "\t" << u128_to_hex(b_leaves[i]) << "\n";
      
    // std::cout << "\n== clo_root first 10 ==\n";
    // for (int i = 0; i < (N < 10 ? N : 10); i++)
    //     std::cout << i << "\t" << u128_to_hex(clo_root[i]) << "\n";

    // std::cout << "\n== chi_root first 10 ==\n";
    // for (int i = 0; i < (N < 10 ? N : 10); i++)
    //     std::cout << i << "\t" << u128_to_hex(chi_root[i]) << "\n";

    // std::cout << "\n== b_root first 10 ==\n";
    // for (int i = 0; i < (N < 10 ? N : 10); i++)
    //     std::cout << i << "\t" << u128_to_hex(b_root[i]) << "\n";

    // std::cout << "\n== c_root first 10 ==\n";
    // for (int i = 0; i < (N < 10 ? N : 10); i++)
    //     std::cout << i << "\t" << u128_to_hex(c_root[i]) << "\n";

    // bool ok = true;
    // for (int i = 0; i < N; i++) {
    //     if (b_root[i] != c_root[i]) {
    //         std::cout << "[TB] mismatch at i=" << i
    //                   << " b_root=" << u128_to_hex(b_root[i])
    //                   << " c_root=" << u128_to_hex(c_root[i]) << "\n";
    //         ok = false;
    //         break;
    //     }
    // }

    // if (ok) {
    //     std::cout << "\n[TB] PASS: b_root == c_root for all " << N << " entries.\n";
    //     return 0;
    // } else {
    //     std::cout << "\n[TB] FAIL.\n";
    //     return 2;
    // }
}