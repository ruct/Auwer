#include <random>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

using std::vector;
using std::string;
using std::cout;
using std::endl;

std::mt19937 GGEN(3711);
const int WH = 512;
const int SIZE = 64;
const string NAME = "data_test1_blank";
const int NL = 2400, NR = 2699;
const int CNTSQ = WH/SIZE;
typedef double ld;

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

//  reading picture
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

//  computations
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

    inline ld metric(int a[3][SIZE], int b[3][SIZE]) {
        ld result = 0;
        for (int c = 0; c < 3; ++c)
        for (int T = 0; T < SIZE; ++T)
            result += ld(sqr(a[c][T]-b[c][T]));
        return sqrt(result);
    }
    void compute() {
        for (int p1 = 0; p1 < CNTSQ*CNTSQ; ++p1)
        for (int p2 = 0; p2 < CNTSQ*CNTSQ; ++p2) {
            if (p1 == p2) continue;

            int i1 = p1/CNTSQ, j1 = p1%CNTSQ;
            int i2 = p2/CNTSQ, j2 = p2%CNTSQ;

            ld cur_diss;
            int a[3][SIZE], b[3][SIZE];

            for (int c = 0; c < 3; ++c)
            for (int i = 0; i < SIZE; ++i) {
                a[c][i] = PICTURE[i1][j1].p[c][i][SIZE-1];
                b[c][i] = PICTURE[i2][j2].p[c][i][0];
            }
            cur_diss = metric(a, b);
            diss[0][p1][p2] = diss[1][p2][p1] = cur_diss;

            for (int c = 0; c < 3; ++c)
            for (int j = 0; j < SIZE; ++j) {
                a[c][j] = PICTURE[i1][j1].p[c][SIZE-1][j];
                b[c][j] = PICTURE[i2][j2].p[c][0][j];
            }
            cur_diss = metric(a, b);
            diss[2][p1][p2] = diss[3][p2][p1] = cur_diss;
        }
    }
}

//  chromosome
namespace ga {
    struct chromo {
        int perm[CNTSQ][CNTSQ];
        int wher[CNTSQ*CNTSQ];
        mutable ld eval = -1;

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
                out << endl;
            }
        }

        ld fit() const {
            if (eval != -1)
                return eval;

            eval = 0;
            for (int i = 0; i < CNTSQ; ++i)
            for (int j = 0; j < CNTSQ; ++j) {
                if (j+1 < CNTSQ)
                    eval += diss[0][perm[i][j]][perm[i][j+1]];
                if (i+1 < CNTSQ)
                    eval += diss[2][perm[i][j]][perm[i+1][j]];
            }
        }
    };

    void random_chromo(chromo& trg) {
        vector <int> cont(CNTSQ*CNTSQ);
        std::iota(cont.begin(), cont.end(), 0);
        std::shuffle(cont.begin(), cont.end(), GGEN);

        for (int i = 0; i < CNTSQ; ++i)
        for (int j = 0; j < CNTSQ; ++j)
            trg.perm[i][j] = cont[i*CNTSQ+j];
        trg.make_wher();
    }
}


const int CORES = 4;
const int GENS = 100;
const int POP = 1000;
const int PREAP = 40;

//  population & kernel
namespace ga {
    struct population {
        int nogen;
        vector <chromo> stock;

        void init() {
            nogen = 0;
            stock.resize(POP);

            for (int i = 0; i < POP; ++i)
                random_chromo(stock[i]);
        }
    };
    struct kernel {
        static const int TCNTSQ = 2*CNTSQ+1;

    private:
        int located = 0;
        vector <int> edge, plan;
        int perm[TCNTSQ][TCNTSQ];
        int mni = TCNTSQ, mxi = -1;
        int mnj = TCNTSQ, mxj = -1;

    public:
        void write();
        kernel();

        //  assuming i in [0, TCNTSQ), j in [0, TCNTSQ)
        int get(int i, int j) const { return perm[i][j]; }
        bool can(int i, int j) const {
            if (perm[i][j] != -1)
                return false;

            if (mxi-mni+1 == CNTSQ && (i < mni || i > mxi))
                return false;
            if (mxj-mnj+1 == CNTSQ && (j < mnj || j > mxj))
                return false;
            return true;
        }
        void place(int i, int j) {

        }
    };
    kernel::kernel() {
        for (int i = 0; i < TCNTSQ; ++i)
        for (int j = 0; j < TCNTSQ; ++j)
            perm[i][j] = -1;
    }
    void kernel::write() {
        for (int i = 0; i < TCNTSQ; ++i) {
            for (int j = 0; j < TCNTSQ; ++j)
                cout << perm[i][j] << "\t";
            cout << endl;
        }
    }

}

//  genetic operators
namespace ga {
    void mutate(chromo& t, std::mt19937& gen);
    void cross(chromo const& a, chromo const& b, std::mt19937& gen);


}
int main(){

}
