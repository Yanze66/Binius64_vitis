#include <iostream>
#include <iomanip>
#include <cstdint>

using u8 = uint8_t;

// 
u8 clmul4_comb(u8 a, u8 b) {
    u8 acc = 0;
    for (int i = 0; i < 4; i++) {
        if ((b >> i) & 1) {
            acc ^= (a << i);
        }
    }
    return acc;
}

int main() {
    std::cout << "static const u8 table[256] = {\n";

    for (int a = 0; a < 16; a++) {
        std::cout << "    // a = 0x" 
                  << std::hex << std::setw(1) << a << std::dec << "\n    ";

        for (int b = 0; b < 16; b++) {
            u8 val = clmul4_comb(a, b);

            std::cout << "0x"
                      << std::hex << std::setw(2) << std::setfill('0')
                      << (int)val << std::dec;

            if (!(a == 15 && b == 15)) std::cout << ",";
        }

        std::cout << "\n";
    }

    std::cout << "};\n";

    return 0;
}
