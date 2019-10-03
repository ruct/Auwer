#include <random>
#include <fstream>
#include <iostream>

using std::string;

std::mt19937 gen(3711);
const int WH = 512;
const int size = 64;
const int NL = 1200, NR = 2700;
const string name = "data_train";
const int cntsq = WH/size;
int read_int(std::ifstream& in) {
    int val;
    in.read(reinterpret_cast <char*>(&val), sizeof(val));
    return val;
}
void write_int(std::ofstream& out, int val) {
    out.write(reinterpret_cast <char*>(&val), sizeof(val));
}

struct square {
    uint8_t r[size][size],
        g[size][size],
        b[size][size];
    square() {
        for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            r[i][j] = g[i][j] = b[i][j] = 0;
    }
    square(square const& o) {
        for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j) {
            r[i][j] = o.r[i][j];
            g[i][j] = o.g[i][j];
            b[i][j] = o.b[i][j];
        }
    }
    void read(std::ifstream& in) {
        for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j) {
            int cur = read_int(in);
            r[i][j] = (cur>>16),
            g[i][j] = ((cur>>8)&255),
            b[i][j] = (cur&255);
        }
    }
    void write(std::ofstream& out) {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j)
                std::cout << r[i][j] << " " << g[i][j] << " " << b[i][j] << "\n";
            std::cout << "\n";
        }
    }
};
typedef double ld;
#define sqr(a) ((a)*(a))

ld horisubs(square const& lft, square const& rht) {
    ld f[] = {0, 0, 0},
        s[] = {0, 0, 0},
        t[] = {0, 0, 0};
    for (int i = 0; i < size; ++i) {
        ld clft[] = {lft.r[i][size-1], lft.g[i][size-1], lft.b[i][size-1]},
            crht[] = {rht.r[i][0], rht.g[i][0], rht.b[i][0]};

        for (int j = 0; j < 3; ++j)
            f[j] += clft[j]*crht[j];
        for (int j = 0; j < 3; ++j)
            s[j] += sqr(clft[j]);
        for (int j = 0; j < 3; ++j)
            t[j] += sqr(crht[j]);
    }
    for (int i = 0; i < 3; ++i)
        s[i] = sqrtl(s[i]);
    for (int i = 0; i < 3; ++i)
        t[i] = sqrtl(t[i]);


    ld res[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        res[i] = f[i]/(s[i]*t[i]);

    std::cout << res[0] << " " << res[1] << " " << res[2] << std::endl;
    ld sim = sqr(res[0])+sqr(res[1])+sqr(res[2]);
    return sim;
}

square picture[cntsq][cntsq];
void read_picture(std::ifstream& in) {
    for (int i = 0; i < cntsq; ++i)
    for (int j = 0; j < cntsq; ++j)
        picture[i][j].read(in);
}
int main() {
    std::ios::sync_with_stdio(false);
    string binput = R"(C:\Users\Main\Base\huawei\parsed\)"+name+"\\"+std::to_string(size);

    std::cout.precision(6);
    for (int n = NL; n < NR; ++n) {
        std::cout << "picture #" << n << std::endl;
        string fbin = binput+"\\"+std::to_string(n)+".txt";
        std::ifstream img_in(fbin, std::ios::binary);

        read_picture(img_in);
        while (true) {
            char type;
            std::cin >> type;
            if (type == 'q')
                break;

            int i1, j1, i2, j2;
            std::cin >> i1 >> j1 >> i2 >> j2;
            std::cout << std::fixed << horisubs(picture[i1][j1], picture[i2][j2]) << std::endl;
        }
    }
}
