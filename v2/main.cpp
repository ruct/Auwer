#include <random>
#include <fstream>
#include <assert.h>
#include <iostream>

using std::string;

std::mt19937 gen(3711);
const int WH = 512;
const int size = 64;
const int NL = 1345, NR = 2700;
const string NAME = "data_train";
const int cntsq = WH/size;

namespace augm {
    int read_int(std::ifstream&);
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
            int cur = augm::read_int(in);
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


square picture[cntsq][cntsq];
namespace augm {
    int read_int(std::ifstream& in) {
        int val;
        in.read(reinterpret_cast <char*>(&val), sizeof(val));
        return val;
    }
    void write_int(std::ofstream& out, int val) {
        out.write(reinterpret_cast <char*>(&val), sizeof(val));
    }
    void read_picture(std::ifstream& in) {
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            picture[i][j].read(in);
    }
}

typedef double ld;
#define sqr(a) ((a)*(a))


int16_t dfans_perm[cntsq][cntsq];
namespace augm {
    bool read_default_answer(std::ifstream& in, int look = -1) {
        string look_for = std::to_string(look);
        while (look_for.size() < 4)
            look_for = "0"+look_for;
        if (look == -1)
            look_for.clear();
        look_for += ".png";

        string cur_line;
        bool found = false;
        while (std::getline(in, cur_line))
            if (cur_line.find(look_for) != string::npos) {
                found = true;
                break;
            }
        if (!found)
            return false;

        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            in >> dfans_perm[i][j];

        return true;
    }
    void write_default_answer(std::ostream& trg) {
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            trg << dfans_perm[i][j] << " ";
        trg << "\n";
    }
}

struct chromo {
    int16_t perm[cntsq][cntsq];
    int16_t wher[cntsq*cntsq];
    ld eval = -1;

    void mk_wher() {
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            wher[perm[i][j]] = i*cntsq+j;
    }
    void shftL() {
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j+1 < cntsq; ++j) {
            -- wher[perm[i][j+1]];
            perm[i][j] = perm[i][j+1];
        }
    }
    void shftR() {
        for (int i = 0; i < cntsq; ++i)
        for (int j = cntsq-1; i > 0; --i) {
            ++ wher[perm[i][j-1]];
            perm[i][j] = perm[i][j-1];
        }
    }
    void shftD() {
        for (int j = 0; j < cntsq; ++j)
        for (int i = cntsq-1; i > 0; --i) {
            wher[perm[i-1][j]] += cntsq;
            perm[i][j] = perm[i-1][j];
        }
    }
    void shftU() {
        for (int j = 0; j < cntsq; ++j)
        for (int i = 0; i+1 < cntsq; ++i) {
            wher[perm[i+1][j]] -= cntsq;
            perm[i][j] = perm[i+1][j];
        }
    }

    void shft(int d) {
        if (d == 0) shftR();
        if (d == 1) shftL();
        if (d == 2) shftD();
        if (d == 3) shftU();
    }
};

//  GA_augmented
namespace GA {
    int dx[] = {0, 0, 1, -1},
        dy[] = {1, -1, 0, 0};
    //L->R, R->L, U->D, D->U;

    //  Initial coordinates;
    ld fitLR[cntsq*cntsq][cntsq*cntsq];
    ld fitUD[cntsq*cntsq][cntsq*cntsq];

    void clear() {
        for (int i = 0; i < cntsq*cntsq; ++i)
        for (int j = 0; j < cntsq*cntsq; ++j) {
            fitLR[i][j] = -1;
            fitUD[i][j] = -1;
        }
    }

    ld horisubs(square const& lft, square const& rht) {
        ld diff[] = {0, 0, 0};
        for (int i = 0; i < size; ++i) {
            diff[0] += sqr(ld(lft.r[i][size-1])-ld(rht.r[i][0]));
            diff[1] += sqr(ld(lft.g[i][size-1])-ld(rht.g[i][0]));
            diff[2] += sqr(ld(lft.b[i][size-1])-ld(rht.b[i][0]));
        }

        ld res = sqrtl(diff[0]+diff[1]+diff[2]);
        return res;
    }
    ld vertsubs(square const& up, square const& dw) {
        ld diff[] = {0, 0, 0};
        for (int j = 0; j < size; ++j) {
            diff[0] += sqr(ld(up.r[size-1][j]-ld(dw.r[0][j])));
            diff[1] += sqr(ld(up.g[size-1][j]-ld(dw.g[0][j])));
            diff[2] += sqr(ld(up.b[size-1][j]-ld(dw.b[0][j])));
        }

        ld res = sqrtl(diff[0]+diff[1]+diff[2]);
        return res;
    }

    ld pair_match(square const& a, square const& b, int dir) {
        if (dir == 0) return horisubs(a, b);
        else if (dir == 1) return horisubs(b, a);
        else if (dir == 2) return vertsubs(a, b);
        else return vertsubs(b, a);
    }

    void recalc() {
        for (int fst = 0; fst < cntsq*cntsq; ++fst)
        for (int snd = 0; snd < cntsq*cntsq; ++snd) {
            if (fst == snd) continue;
            fitLR[fst][snd] = horisubs(picture[fst/cntsq][fst%cntsq],
                                       picture[snd/cntsq][snd%cntsq]);
            fitUD[fst][snd] = vertsubs(picture[fst/cntsq][fst%cntsq],
                                       picture[snd/cntsq][snd%cntsq]);
        }
    }
    //  Initial coordinates;
    inline ld pair_match(int fst, int snd, int dir) {
        if (dir&1) std::swap(fst, snd);
        if (dir < 0) return fitLR[fst][snd];
        else return fitUD[fst][snd];
    }

    ld fit_cost(chromo& o) {
        if (o.eval != -1)
            return o.eval;

        ld eval = 0;
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j) {
            if (j+1 < cntsq)
                eval += fitLR[o.perm[i][j]][o.perm[i][j+1]];
            if (i+1 < cntsq)
                eval += fitUD[o.perm[i][j]][o.perm[i+1][j]];
        }
        return o.eval = eval;
    }
}

#include <algorithm>
//  Ga_algorithm
const int FILTER = 40;
chromo scent[FILTER];
namespace GA {
    void stupid_shuffle(chromo& trg) {
        std::vector <uint16_t> lol(cntsq*cntsq);
        std::iota(lol.begin(), lol.end(), 0);
        std::shuffle(lol.begin(), lol.end(), gen);

        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            trg.perm[i][j] = lol[i*cntsq+j];
        trg.mk_wher();
        trg.eval = -1;
    }
    void few_swaps(chromo& trg, int swaps) {
        for (int T = 0; T < swaps; ++T) {
            int i1 = gen()%cntsq, j1 = gen()%cntsq,
                i2 = gen()%cntsq, j2 = gen()%cntsq;

            std::swap(trg.wher[trg.perm[i1][j1]],
                    trg.wher[trg.perm[i2][j2]]);
            std::swap(trg.perm[i1][j1], trg.perm[i2][j2]);
        }
    }

//    void init() {
//        for (int i = 0; i < FILTER; ++i) {
//            stupid_shuffle(scent[i]);
//        }
//    }
}

int main() {
    std::ios::sync_with_stdio(false);
    string shimg_prefix = R"(C:\Users\Main\Base\huawei\parsed\)"+NAME+"\\"+std::to_string(size);
    string crimg_prefix = R"(C:\Users\Main\Base\huawei\)"+NAME+std::to_string(size);
    string dfans_prefix = R"(C:\Users\Main\Base\huawei\)"+NAME+"\\"+NAME+"_"+std::to_string(size)+"_answers.txt";

    std::cout.precision(6);
    for (int n = NL; n < NR; ++n) {
        std::cout << "picture #" << n << std::endl;
        string shimg = shimg_prefix+"\\"+std::to_string(n)+".txt";
        string crimg = crimg_prefix+"\\"+std::to_string(n)+".txt";

        std::ifstream shimg_in(shimg, std::ios::binary);
        std::ifstream crimg_in(crimg);

        augm::read_picture(shimg_in);

///aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa///

        std::cout << dfans_prefix << "\n";
        std::ifstream dfans_in(dfans_prefix);
        while (true) {
            int look = n;

            std::cout << dfans_prefix << " " << n << std::endl;
            system("pause");

            bool success = augm::read_default_answer(dfans_in, look);
            std::cout << std::boolalpha << success << "\n";

        }
    }
}
