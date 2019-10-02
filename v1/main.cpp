#include <string>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <algorithm>

typedef long long ll;
#define sqr(a) ((a)*(a))
int read_int(std::ifstream& in) {
    int val;
    in.read(reinterpret_cast<char *>(&val), sizeof(val));
    return val;
}
void write_int(std::ofstream& out, int val) {
    out.write(reinterpret_cast<char*>(&val), sizeof(val));
}

std::mt19937 gen(3711);
const int WH = 512;
const int size = 64;
const int NL = 1200, NR = 2700;
const string name = "data_train";
const int cntsq = WH/size;

ll best_estim = 9e18;
const int Tries = 2e5;
uint16_t P[cntsq][cntsq];
uint16_t BP[cntsq][cntsq];
uint8_t r[WH][WH], g[WH][WH], b[WH][WH];
void read(std::ifstream& in) {
    for (int i = 0; i < WH; ++i)
    for (int j = 0; j < WH; ++j) {
        int mask = read_int(in);
//        std::cout << mask << std::endl;
//        system("pause");

        r[i][j] = mask>>16,
        g[i][j] = (mask>>8)&255,
        b[i][j] = mask&255;
//        std::cout << r[i][j] << " " << g[i][j] << " " << b[i][j] << std::endl;
    }
}

int dist(int i1, int j1, int i2, int j2) {
    if (std::min({i1, i2, j1, j2}) < 0 || std::max({i1, i2, j1, j2}) >= WH) {
            std::cout << i1 << " " << j1 << " " << i2 << " " << j2 << std::endl;
        std::cout << "haahahhahah" << std::endl;
        system("pause");
    }

    return sqr(r[i1][j1]-r[i2][j2])
        +sqr(g[i1][j1]-g[i2][j2])
        +sqr(b[i1][j1]-b[i2][j2]);
}
ll estim(uint16_t perm[cntsq][cntsq]) {
    ll res = 0;
//    std::cout << "Estimating" << std::endl;
//    for (int i = 0; i < cntsq; ++i) {
//        for (int j = 0; j < cntsq; ++j)
//            std::cout << perm[i][j] << " ";
//        std::cout << "\n";
//    }
//    std::cout << std::endl << std::endl;
//    system("pause");


    for (int sqi = 0; sqi < cntsq; ++sqi)
    for (int sqj = 0; sqj < cntsq; ++sqj) {
//        std::cout << sqi << " " << sqj << " " << res << std::endl;
        int f = perm[sqi][sqj];
        int oqi1 = f/cntsq, oqj1 = f%cntsq;

//        std::cout << "calced" << std::endl;
//        std::cout << oqi1 << " " << oqj1 << " " << oqi2 << " " << oqj2 << std::endl;
        int i1 = oqi1*size,
            j1 = oqj1*size;

        if (sqi+1 < cntsq) {
            int s = perm[sqi+1][sqj];
            int oqi2 = s/cntsq, oqj2 = s%cntsq;

            int i2 = oqi2*size,
                j2 = oqj2*size;

//                std::cout << "LOOOL " << i1+size-1 << " " << j1 << " " << i2 << " " << j2 << std::endl;
            for (int j = 0; j < size; ++j)
                res += dist(i1+size-1, j1+j, i2, j2+j);
        }
        if (sqj+1 < cntsq) {
            int s = perm[sqi][sqj+1];
            int oqi2 = s/cntsq, oqj2 = s%cntsq;


            int i2 = oqi2*size,
                j2 = oqj2*size;

//                std::cout << "LOOOLKKEEEk " << i1 << " " << j1+size-1 << " " << i2 << " " << j2 << std::endl;
            for (int i = 0; i < size; ++i)
                res += dist(i1+i, j1+size-1, i2+i, j2);
        }
    }
//    std::cout << "mda lol " << res << std::endl;
    return res;
}

void stupid_shuffle(uint16_t perm[cntsq][cntsq]) {
    std::vector <uint16_t> lol;
    lol.resize(cntsq*cntsq);
    std::iota(lol.begin(), lol.end(), 0);
    std::shuffle(lol.begin(), lol.end(), gen);

//    std::cout << "AHAHAHDIAUHDAHIAHDAHID " << std::endl;
    for (int i = 0; i < cntsq; ++i)
    for (int j = 0; j < cntsq; ++j)
        perm[i][j] = lol[i*cntsq+j];
}
void few_swaps(uint16_t perm[cntsq][cntsq], int swaps) {
    for (int T = 0; T < swaps; ++T) {
        int i1 = gen()%cntsq, j1 = gen()%cntsq,
            i2 = gen()%cntsq, j2 = gen()%cntsq;
        std::swap(perm[i1][j1], perm[i2][j2]);
    }
}

void copy(uint16_t Psrc[cntsq][cntsq], uint16_t Ptrg[cntsq][cntsq]) {
    for (int i = 0; i < cntsq; ++i)
    for (int j = 0; j < cntsq; ++j)
        Ptrg[i][j] = Psrc[i][j];
}

void solve() {
    for (int i = 0; i < cntsq; ++i)
    for (int j = 0; j < cntsq; ++j)
        P[i][j] = i*cntsq+j;
//    std::cout << "kek" << std::endl;

    best_estim = estim(P);
    stupid_shuffle(P);

//    std::cout << "Was " << best_estim << "\n";
    copy(P, BP);

    for (int T = 0; T < Tries; ++T) {
        if (__builtin_popcount(T) == 1)
            std::cout << T << " " << best_estim << std::endl;
        copy(P, BP);
        few_swaps(BP, 4);
//        stupid_shuffle(BP);

        ll cur_estim = estim(BP);
        if (cur_estim < best_estim) {
            best_estim = cur_estim;
            copy(BP, P);
        }
//        std::cout << T << std::endl;
//        if (T%50 == 0) {
//            std::cout << "CUR " << best_estim << std::endl;
//            system("pause");
//        }
    }
//    std::cout << "Became " << best_estim << "\n";
}

std::string get_number(int n) {
    std::string res = std::to_string(n);
    while (res.size() < 4)
        res = "0"+res;
    return res;
}

const int IM = 3000;
ll best_in_file[IM];
std::vector <int> best_perms[IM];
void write(std::ofstream& out, int cur_id) {
    if (best_estim >= )
    out << get_number(cur_id) << "";
    for (int i = 0; i < cntsq; ++i)
    for (int j = 0; j < cntsq; ++j)
        out << P[i][j] << " ";
    out << "\n";
}

//format:
//image_id best_estimation_found
//perm1 perm2 ... perm[cntsq*cntsq]
void make_bests(std::ifstream& in) {
    int im_id;
    while (in >> im_id) {
        ll best_found;
        in >> best_found;

        best_in_file[im_id] = best_found;
        best_perms[im_id].resize(cntsq*cntsq);
        for (int i = 0; i < ctnsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            in >> best_perms[im_id][i*cntsq+j];
    }
}

int main(){
    std::ios::sync_with_stdio(0);
    std::string binput = R"(C:\Users\Main\Base\huawei\parsed\)"+name+"\\"+std::to_string(size);
    std::string output = binput+"\\answers.txt";

    std::ofstream out(output);
    for (int n = NL; n < NL+1; ++n) {
        std::string fbin = binput+"\\"+std::to_string(n)+".txt";
        std::ifstream in(fbin, std::ios::binary);

        read(in);
        solve();
        write(out, n);
    }
    std::cout << best_estim << std::endl;
}
