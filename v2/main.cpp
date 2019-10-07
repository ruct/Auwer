#include <tuple>
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

//  GA_fitness
namespace GA {
    int dx[] = {0, 0, 1, -1},
        dy[] = {1, -1, 0, 0};
    //L->R, R->L, U->D, D->U;

    //  Initial coordinates;
    int bbuddy[4*cntsq*cntsq];
    ld fitLR[cntsq*cntsq][cntsq*cntsq];
    ld fitUD[cntsq*cntsq][cntsq*cntsq];

    void clear() {
        for (int i = 0; i < cntsq*cntsq; ++i) {
            for (int d = 0; d < 4; ++d)
                bbuddy[i*4+d] = 0;
            for (int j = 0; j < cntsq*cntsq; ++j) {
                fitLR[i][j] = -1;
                fitUD[i][j] = -1;
            }
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

    //  Initial coordinates;
    inline ld pair_match(int fst, int snd, int dir) {
        if (dir&1) std::swap(fst, snd);
        if (dir < 0) return fitLR[fst][snd];
        else return fitUD[fst][snd];
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
        for (int fst = 0; fst < cntsq*cntsq; ++fst) {
            for (int d = 0; d < 4; ++d) {
                std::pair <ld, int> mn = {1e18, -1};
                for (int snd = 0; snd < cntsq*cntsq; ++snd)
                    mn = std::min(mn, {pair_match(fst, snd, d), snd});
                bbuddy[fst*4+d] = mn.second;
            }
        }
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
//  Ga_random
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
        trg.eval = -1;
    }

}

const int FILTER = 40;
chromo scent[FILTER];
const int ST_FEW_SWAPS = 50;
const int FSKIP = 10;
namespace GA {

    void init() {
        for (int i = 0; i < FILTER; ++i) {
            stupid_shuffle(scent[i]);
            for (int j = 0; j < ST_FEW_SWAPS; ++j) {
                chromo cur = scent[i];
                few_swaps(cur, 3);
                if (fit_cost(cur) > fit_cost(scent[i]))
                    scent[i] = cur;
            }
        }
    }

    const int tcntsq = 2*cntsq+1;
    struct temp_chromo {
        int mni, mxi, mnj, mxj;
        int16_t perm[tcntsq][tcntsq];
        temp_chromo() {
            for (int i = 0; i < tcntsq; ++i)
            for (int j = 0; j < tcntsq; ++j)
                perm[i][j] = -1;

            mni = mnj = tcntsq;
            mxi = mxj = -1;
        }

        void set(int i, int j, int val) {
            perm[i][j] = val;
            mni = std::min(mni, i),
            mxi = std::max(mxi, i);
            mnj = std::min(mnj, j),
            mxj = std::max(mxj, j);
        }

        bool can(int i, int j) {
            if (perm[i][j] == -1)
                return false;

            int _mni = std::min(mni, i),
                _mxi = std::max(mxi, i);
            int _mnj = std::min(mnj, j),
                _mxj = std::max(mxj, j);

            if (_mxi-_mni+1 > cntsq)
                return false;
            if (_mxj-_mnj+1 > cntsq)
                return false;
            return true;
        }
    };

    int del_rnd(std::vector <int>& v) {
        int i = gen()%v.size();
        std::swap(v[i], v.back());
        i = v.back();
        v.pop_back();
        return i;
    }

    void cross(chromo& a, chromo& b, chromo& t) {
        temp_chromo temp;
        std::vector <int> lft(cntsq*cntsq);
        std::iota(lft.begin(), lft.end(), 0);

        int placed = 1;
        char used[cntsq*cntsq];
        for (int i = 0; i < cntsq*cntsq; ++i)
            used[i] = 0;
        std::vector <std::pair <int, int> > bound = {{cntsq, cntsq}};
        used[temp.perm[cntsq][cntsq] = del_rnd(lft)] = 1;

        auto bound_update = [&dx, &dy, &bound]() {

        };

        while (true) {
            //{index in number, direction, 2-agreed id}
            std::vector <std::tuple <int, int, int> > suitable;

            for (int TT = 0; TT < bound.size(); ++TT) {
                int i = bound[TT].first, j = bound[TT].second;

                int cur_id = temp.perm[i][j];
                int pos1 = a.wher[cur_id],
                    pos2 = b.wher[cur_id];

                for (int d = 0; d < 4; ++d) {
                    if (!temp.can(i+dx[d], j+dy[d]))
                        continue;

                    int ni1 = (pos1/cntsq)+dx[d], nj1 = (pos1%cntsq)+dy[d];
                    int ni2 = (pos2/cntsq)+dx[d], nj2 = (pos2%cntsq)+dy[d];

                    if (ni1 < 0 || ni1 >= cntsq) continue;
                    if (nj1 < 0 || nj1 >= cntsq) continue;
                    if (ni2 < 0 || ni2 >= cntsq) continue;
                    if (nj2 < 0 || nj2 >= cntsq) continue;
                    if (a.perm[ni1][nj1] == b.perm[ni2][nj2] && !used[a.perm[ni1][nj1]])
                        suitable.emplace_back(TT, d, a.perm[ni1][nj1]);
                }
            }

            if (suitable.size() && gen()%100 > FSKIP) {
                int srnd = gen()%suitable.size();
                int TT = std::get<0> (suitable[srnd]);
                int id = std::get<2> (suitable[srnd]);
                int d = std::get<1> (suitable[srnd]);

                int i = bound[TT].first,
                    j = bound[TT].second;

                ++placed;
                used[id] = 1;
                temp.set(i+dx[d], j+dy[d], id);
                lft.erase(std::find(lft.begin(), lft.end(), id));
                bound_update();
                continue;
            }
        }
    }

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
