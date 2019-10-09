#include <random>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

std::mt19937 GEN(3711);
const int WH = 512;
const int SIZE = 64;
const string NAME = "data_test1_blank";
const int NL = 2400, NR = 2699;
const int CNTSQ = WH/SIZE;


namespace augm {
    int read_int(std::ifstream&);
}
struct square {
    int p[3][SIZE][SIZE];
    square() {
        for (int c = 0; c < 3; ++c)
        for (int i = 0; i < SIZE; ++i)
        for (int j = 0; j < SIZE; ++j)
            p[c][i][j] = 0;
    }
    square(square const& o) {
        for (int c = 0; c < 3; ++c)
        for (int i = 0; i < SIZE; ++i)
        for (int j = 0; j < SIZE; ++j)
            p[c][i][j] = o.p[c][i][j];
    }
};

square PICTURE[CNTSQ][CNTSQ];
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
        int picture_buf[3][WH][WH];

        for (int j = 0; j < WH; ++j)
        for (int i = 0; i < WH; ++i) {
            int cur = read_int(in);
            int bands[] = {cur>>16, (cur>>8)&255, cur&255};
            for (int c = 0; c < 3; ++c)
                picture_buf[c][i][j] = bands[c];
        }
        for (int c = 0; c < 3; ++c)
        for (int i = 0; i < CNTSQ; ++i)
        for (int j = 0; j < CNTSQ; ++j)
            for (int x = 0; x < SIZE; ++x)
            for (int y = 0; y < SIZE; ++y)
                PICTURE[i][j].p[c][x][y] = picture_buf[c][i*SIZE+x][j*SIZE+y];
    }
}



typedef double ld;
namespace ga {
    int dx[] = {0, 0, 1, -1},
        dy[] = {1, -1, 0, 0};
    //L->R, R->L, U->D, D->U;
    int bbuddy[4][CNTSQ*CNTSQ];
    ld diss[4][CNTSQ*CNTSQ][CNTSQ*CNTSQ];

    void reset() {
        for (int d = 0; d < 4; ++d)
        for (int i = 0; i < CNTSQ*CNTSQ; ++i) {
            bbuddy[d][i] = -1;
            for (int j = 0; j < CNTSQ*CNTSQ; ++j)
                diss[d][i][j] = -1;
        }
    }

#define sqr(a) ((a)*(a))
#define cbr(a) ((a)*(a)*(a))

    inline ld metric(std::vector <int> const& a,
                    std::vector <int> const& b) {
        ld result = 0;
        for (int T = 0; T < a.size(); ++T)
            result += ld(sqr(a[T]-b[T]));
        return sqrt(result);
    }

}

namespace ga {
    struct chromo {
        int perm[CNTSQ][CNTSQ];
        int wher[CNTSQ*CNTSQ];
        ld eval = -1;

        bool operator!=(chromo const& o) const {
            for (int i = 0; i < CNTSQ; ++i)
            for (int j = 0; j < CNTSQ; ++j)
                if (perm[i][j] != o.perm[i][j])
                    return true;
            return false;
        }
        void make_wher() {
            for (int i = 0; i < CNTSQ; ++i)
            for (int j = 0; j < CNTSQ; ++j)
                wher[perm[i][j]] = i*CNTSQ+j;
        }
        void write(std::ofstream& out) const {
            for (int i = 0; i < CNTSQ; ++i) {
                for (int j = 0; j < CNTSQ; ++j)
                    out << perm[i][j] << "\t";
                out << std::endl;
            }
        }

    };

}

int main(){

}
