#include <tuple>
#include <time.h>
#include <chrono>
#include <future>
#include <random>
#include <fstream>
#include <assert.h>
#include <iostream>
#include <windows.h>

using std::string;

std::mt19937 gen(3711);
const int WH = 512;
const int size = 64;
const int NL = 1200, NR = 1220;
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

int16_t RPICTURE[WH][WH];
int16_t GPICTURE[WH][WH];
int16_t BPICTURE[WH][WH];
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
        for (int j = 0; j < WH; ++j)
        for (int i = 0; i < WH; ++i) {
            int cur = read_int(in);
            int r = cur>>16,
                g = (cur>>8)&255,
                b = cur&255;
            RPICTURE[i][j] = r;
            GPICTURE[i][j] = g;
            BPICTURE[i][j] = b;
        }
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            for (int x = 0; x < size; ++x)
            for (int y = 0; y < size; ++y) {
                picture[i][j].r[x][y] = RPICTURE[i*size+x][j*size+y];
                picture[i][j].g[x][y] = GPICTURE[i*size+x][j*size+y];
                picture[i][j].b[x][y] = BPICTURE[i*size+x][j*size+y];
            }
    }
}

typedef double ld;
#define sqr(a) ((a)*(a))
#define cbr(a) ((a)*(a)*(a))


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

    bool operator!=(chromo const& o) const {
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            if (perm[i][j] != o.perm[i][j])
                return true;
        return false;
    }
    void mk_wher() {
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            wher[perm[i][j]] = i*cntsq+j;
    }
    void write() const {
        for (int i = 0; i < cntsq; ++i) {
            for (int j = 0; j < cntsq; ++j)
                std::cout << perm[i][j] << "\t";
            std::cout << std::endl;
        }
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
                bbuddy[i*4+d] = -1;
            for (int j = 0; j < cntsq*cntsq; ++j) {
                fitLR[i][j] = -1;
                fitUD[i][j] = -1;
            }
        }
    }

    ld horisubs(square const& lft, square const& rht) {
        ld diff[] = {0, 0, 0};
        for (int i = 0; i < size; ++i) {
            diff[0] += cbr(abs(ld(lft.r[i][size-1])-ld(rht.r[i][0])));
            diff[1] += cbr(abs(ld(lft.g[i][size-1])-ld(rht.g[i][0])));
            diff[2] += cbr(abs(ld(lft.b[i][size-1])-ld(rht.b[i][0])));
//            diff[0] += sqr(ld(lft.r[i][size-1])-ld(rht.r[i][0]));
//            diff[1] += sqr(ld(lft.g[i][size-1])-ld(rht.g[i][0]));
//            diff[2] += sqr(ld(lft.b[i][size-1])-ld(rht.b[i][0]));
        }

        ld res = cbrt(diff[0]+diff[1]+diff[2]);
//        std::cout << "REES1 " << res << " " << diff[0] << " " << diff[1] << " " << diff[2] << std::endl;
        assert(res >= 0);
//        ld res = sqrtl(diff[0]+diff[1]+diff[2]);
        return res;
    }
    ld vertsubs(square const& up, square const& dw) {
        ld diff[] = {0, 0, 0};
        for (int j = 0; j < size; ++j) {
            diff[0] += cbr(abs(ld(up.r[size-1][j])-ld(dw.r[0][j])));
            diff[1] += cbr(abs(ld(up.g[size-1][j])-ld(dw.g[0][j])));
            diff[2] += cbr(abs(ld(up.b[size-1][j])-ld(dw.b[0][j])));
//            diff[0] += sqr(ld(up.r[size-1][j]-ld(dw.r[0][j])));
//            diff[1] += sqr(ld(up.g[size-1][j]-ld(dw.g[0][j])));
//            diff[2] += sqr(ld(up.b[size-1][j]-ld(dw.b[0][j])));
        }

        ld res = cbrt(diff[0]+diff[1]+diff[2]);
//        std::cout << "REES2 " << res << " " << diff[0] << " " << diff[1] << " " << diff[2] << std::endl;
        assert(res >= 0);
//        ld res = sqrtl(diff[0]+diff[1]+diff[2]);
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
        GA::clear();
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
                    if (fst != snd)
                        mn = std::min(mn, {pair_match(fst, snd, d), snd});
                bbuddy[4*fst+d] = mn.second;
            }
        }
    }


    ld fit_cost(chromo& o) {
        if (o.eval != -1)
            return o.eval;

        ld eval = 0;
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j) {
            if (j+1 < cntsq) {
                eval += fitLR[o.perm[i][j]][o.perm[i][j+1]];
            }
            if (i+1 < cntsq) {
                eval += fitUD[o.perm[i][j]][o.perm[i+1][j]];
            }
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

const int CORES = 4;
const int GENS = 100;
const int REAP = 10;
const int PEACE = 2000;
chromo scent[PEACE];
const int ST_FEW_SWAPS = 3;
const int FSKIP = 6;
const int TSKIP = 6;
const int CHANCE_HARE = 0;
const int MUTATION_RATE = 0;
const int FIRST_PHASE_AFTER = 10;

struct TimeLogger {
    ld start;
    ld* target;
    TimeLogger() = delete;
    explicit TimeLogger(ld* target):
        target(target), start(clock()) {};
    ~TimeLogger() {
        *target += ld(clock()-start)/ld(CLOCKS_PER_SEC);
    }
};

int CUR_GEN;
chromo GLOBAL;
namespace GA {

    void init() {
        for (int i = 0; i < PEACE; ++i) {
            stupid_shuffle(scent[i]);
            for (int j = 0; j < ST_FEW_SWAPS; ++j) {
                chromo cur = scent[i];
                few_swaps(cur, 2);
                if (fit_cost(cur) > fit_cost(scent[i]))
                    scent[i] = cur;
            }
        }
        GLOBAL = scent[0];
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
            if (perm[i][j] != -1)
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

        void write() {
            std::cout << "Writing temp_chromo:\n";
            for (int i = 0; i < tcntsq; ++i) {
                for (int j = 0; j < tcntsq; ++j)
                    std::cout << perm[i][j] << "\t";
                std::cout << std::endl;
            }
        }

    };

    void mutate(chromo& t) {

    }
    void cross(chromo const& a, chromo const& b, chromo& t, std::mt19937& genr) {
        temp_chromo temp;

        std::vector <int> lft(cntsq*cntsq);
        std::iota(lft.begin(), lft.end(), 0);

        int placed = 0;
        char used[cntsq*cntsq];
        for (int i = 0; i < cntsq*cntsq; ++i)
            used[i] = 0;
        std::vector <std::pair <int, int> > bound;

        auto bound_update = [&temp, &dx, &dy, &bound]() {
            std::vector <std::pair <int, int> > buf;
            for (auto P : bound) {
                int i = P.first,
                    j = P.second;

                int cnt = 0;
                for (int d = 0; d < 4; ++d) {
                    int ni = i+dx[d],
                        nj = j+dy[d];
                    cnt += temp.can(ni, nj);
                }
                if (cnt > 0) {
                    buf.emplace_back(i, j);
                }
            }
            bound = buf;
        };
        auto place_fragment = [&temp, &placed, &used, &lft, &bound, &bound_update]
            (int i, int j, int id) {
                ++placed;
                used[id] = 1;
//                std::cout << "F" << std::endl;
                temp.set(i, j, id);
//                std::cout << "S" << " " << i << " " << j << std::endl;

                bool happened = 0;
                for (auto it = lft.begin(); it != lft.end(); ++it)
                    if (*it == id) {
                        happened = 1;
                        lft.erase(it);
                        break;
                    }
                assert(happened);
                bound.emplace_back(i, j);
                bound_update();
            };

//        std::cout << "A:\n";
//        a.write();
//        std::cout << "\nB:\n";
//        b.write();
//        std::cout << std::endl;

        place_fragment(cntsq, cntsq, genr()%(cntsq*cntsq));
        while (placed < cntsq*cntsq) {
//            std::cout << placed << std::endl;
//            temp.write();
//            system("pause");

            //{index in number, direction, suitable_id}
            std::vector <std::tuple <int, int, int> > suitable;
            ///  first phase
            if (CUR_GEN > FIRST_PHASE_AFTER){

                for (int TT = 0; TT < bound.size(); ++TT) {
                    int i = bound[TT].first, j = bound[TT].second;

                    int cur_id = temp.perm[i][j];
                    int pos1 = a.wher[cur_id],
                        pos2 = b.wher[cur_id];
                    assert(a.perm[pos1/cntsq][pos1%cntsq] == cur_id);
                    assert(b.perm[pos2/cntsq][pos2%cntsq] == cur_id);

                    for (int d = 0; d < 4; ++d) {
                        if (!temp.can(i+dx[d], j+dy[d]))
                            continue;

                        int ni1 = (pos1/cntsq)+dx[d], nj1 = (pos1%cntsq)+dy[d];
                        int ni2 = (pos2/cntsq)+dx[d], nj2 = (pos2%cntsq)+dy[d];

                        if (ni1 < 0 || ni1 >= cntsq) continue;
                        if (nj1 < 0 || nj1 >= cntsq) continue;
                        if (ni2 < 0 || ni2 >= cntsq) continue;
                        if (nj2 < 0 || nj2 >= cntsq) continue;
                        if (a.perm[ni1][nj1] == b.perm[ni2][nj2] && !used[a.perm[ni1][nj1]]) {
                            suitable.emplace_back(TT, d, a.perm[ni1][nj1]);
                        }
                    }
                }

//                std::cout << "First " << suitable.size() << std::endl;
                if (suitable.size()) {
                    if (genr()%100 < FSKIP) {
                        int dirs[] = {0, 1, 2, 3};
                        std::shuffle(dirs, dirs+4, genr);

                        int P = genr()%bound.size();
                        int i = bound[P].first,
                            j = bound[P].second;
                        for (int _d = 0; _d < 4; ++_d) {
                            int d = dirs[_d];
                            int ni = i+dx[d],
                                nj = j+dy[d];
                            if (!temp.can(ni, nj))
                                continue;
                            place_fragment(ni, nj, lft[genr()%lft.size()]);
                        }
                        continue;
                    }

                    int srnd = genr()%suitable.size();
                    int TT = std::get<0> (suitable[srnd]);
                    int id = std::get<2> (suitable[srnd]);
                    int d = std::get<1> (suitable[srnd]);

                    int i = bound[TT].first,
                        j = bound[TT].second;

                    place_fragment(i+dx[d], j+dy[d], id);

                    continue;
                }
            }

            ///second phase
            {
//                TimeLogger second_phase_time(&sphase_time);
                suitable.clear();

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

                        if (std::min(ni1, nj1) >= 0 && std::max(ni1, nj1) < cntsq) {
                            int id1 = a.perm[ni1][nj1];
                            if (!used[id1] && bbuddy[4*cur_id+d] == id1 && bbuddy[4*id1+(d^1)] == cur_id) {
                                suitable.emplace_back(TT, d, id1);
                            }
                        }
                        if (std::min(ni2, nj2) >= 0 && std::max(ni2, nj2) < cntsq) {
                            int id2 = b.perm[ni2][nj2];
                            if (!used[id2] && bbuddy[4*cur_id+d] == id2 && bbuddy[4*id2+(d^1)] == cur_id) {
                                suitable.emplace_back(TT, d, id2);
                            }
                        }
                    }
                }
//                std::cout << "FINISH" << std::endl;

//                std::cout << "Second phase " << suitable.size() << std::endl;
                if (suitable.size()) {
                    int srnd = genr()%suitable.size();
                    int TT = std::get<0> (suitable[srnd]);
                    int id = std::get<2> (suitable[srnd]);
                    int d = std::get<1> (suitable[srnd]);

                    int i = bound[TT].first,
                        j = bound[TT].second;

                    place_fragment(i+dx[d], j+dy[d], id);
                    continue;
                }
            }

//            std::cout << "try third" << std::endl;
            ///third phase
            {
//                TimeLogger third_phase_time(&tphase_time);
                if (genr()%100 < TSKIP) {
                    int dirs[] = {0, 1, 2, 3};
                    std::shuffle(dirs, dirs+4, genr);

                    int P = genr()%bound.size();
                    int i = bound[P].first,
                        j = bound[P].second;
                    for (int _d = 0; _d < 4; ++_d) {
                        int d = dirs[_d];
                        int ni = i+dx[d],
                            nj = j+dy[d];
                        if (!temp.can(ni, nj))
                            continue;
                        place_fragment(ni, nj, lft[genr()%lft.size()]);
                    }
                    continue;
                }

                int TT = genr()%bound.size();
                int i = bound[TT].first,
                    j = bound[TT].second;
                int id = temp.perm[i][j];

                ld mn = 1e9;
                int mn_id = -1, mn_d = -1;
                for (int d = 0; d < 4; ++d) {
                    int ni = i+dx[d],
                        nj = j+dy[d];

                    if (!temp.can(ni, nj))
                        continue;

                    for (int lft_id : lft) {
                        ld cur = pair_match(id, lft_id, d);
                        if (cur < mn) {
                            mn = cur;
                            mn_d = d;
                            mn_id = lft_id;
                        }
                    }
                }
                place_fragment(i+dx[mn_d], j+dy[mn_d], mn_id);
            }
        }
        assert(temp.mxi-temp.mni+1 == cntsq);
        assert(temp.mxj-temp.mnj+1 == cntsq);

        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            t.perm[i][j] = temp.perm[i+temp.mni][j+temp.mnj];

        if (genr()%100 < MUTATION_RATE)
            mutate(t);
        t.eval = -1;
        t.mk_wher();
    }

    chromo parents[REAP];
    void breedLR(int l, int r) {
        std::mt19937 genr(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        for (int i = l; i <= r; ++i) {
//            std::cout << i << std::endl;
            int a = genr()%REAP,
                b = genr()%(REAP-1);
            if (b >= a) ++b;
//            std::cout << a << " " << b << " " << i << std::endl;
            cross(parents[a], parents[b], scent[i], genr);
        }
//        std::cout << "finished " << l << " " << r << std::endl;
    }
    void breed() {
        std::cout << "start breeding" << std::endl;
        for (int i = 0; i < REAP; ++i)
            parents[i] = scent[i];

        int each = PEACE/CORES;
        assert(PEACE%CORES == 0);

        std::future <void> threads[CORES];
        for (int i = 0; i < CORES; ++i)
            threads[i] = std::async(std::launch::async, breedLR, i*each, (i+1)*each-1);
        for (int i = 0; i < CORES; ++i)
            threads[i].get();

        std::cout << "finished breeding" << std::endl;
    }

    void harvest() {
        std::vector <std::pair <ld, int> > _scent(PEACE);
        for (int i = 0; i < PEACE; ++i)
            _scent[i] = {fit_cost(scent[i]), i};
        std::sort(_scent.begin(), _scent.end());

        std::vector <chromo> buf(REAP);

        int sz = 0;
        for (int i = 0; sz < REAP; ++i) {
            if (_scent[i+1].first > _scent[i].first) {
                buf[sz] = scent[_scent[i].second];
                ++sz;
            }
        }
        for (int i = 0; i < REAP; ++i)
            scent[i] = buf[i];

        if (fit_cost(GLOBAL) > fit_cost(scent[0]))
            GLOBAL = scent[0];
        if (gen()%100 < CHANCE_HARE)
            scent[1] = GLOBAL;

        std::cout << "harvested" << std::endl;
        for (int i = 0; i < REAP; ++i)
            std::cout << fit_cost(scent[i]) << " ";
        std::cout << std::endl;
//        std::shuffle(scent, scent+REAP, gen);
    }

    int NOFIMAGE;
    int16_t GA_ans[cntsq][cntsq];
    void genetic_algorithm(std::ofstream& out) {
        init();
        harvest();

        for (CUR_GEN = 0; CUR_GEN < GENS; ++CUR_GEN) {
            breed(), harvest();

            ld sum = 0, mn = 1e18;
            for (int j = 0; j < REAP; ++j) {
                sum += fit_cost(scent[j]);
                mn = std::min(mn, fit_cost(scent[j]));
            }
            std::cout << NOFIMAGE << " Generation #" << CUR_GEN << " " << sum << " " << mn << std::endl;
            std::cout << "MINIMUM " << fit_cost(GLOBAL) << std::endl;
//            system("pause");
        }

        std::cout << "\n\n\n";
//        GLOBAL.write();

        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            GA_ans[i][j] = GLOBAL.perm[i][j];
    }
}


int main() {
    std::ios::sync_with_stdio(false);

    //  shuffled
    string shimg_prefix = R"(C:\Users\Main\Base\huawei\parsed\)"+NAME+"\\"+std::to_string(size);
    //  correct
    string crimg_prefix = R"(C:\Users\Main\Base\huawei\)"+NAME+std::to_string(size);
    //  default answer
    string dfans_prefix = R"(C:\Users\Main\Base\huawei\)"+NAME+"\\"+NAME+"_"+std::to_string(size)+"_answers.txt";

    string ans_write_path = R"(C:\Users\Main\Base\huawei\parsed\)" + NAME + "\\" +
            std::to_string(size) + "\\answers.txt";
    std::ofstream ans_out(ans_write_path);
    for (int n = NL; n <= NL+30; ++n) {
        GA::NOFIMAGE = n;
        string shimg = shimg_prefix+"\\"+std::to_string(n)+".txt";
        string crimg = crimg_prefix+"\\"+std::to_string(n)+".txt";
        std::ifstream shimg_in(shimg, std::ios::binary);
        std::ifstream crimg_in(crimg);

        augm::read_picture(shimg_in);
        GA::recalc();

        GA::genetic_algorithm(ans_out);

        ans_out << std::to_string(n) << ".png" << std::endl;
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            ans_out << GA::GA_ans[i][j] << " ";
        ans_out << std::endl;
        std::cout << "\n";
    }
//    std::cout << dfans_prefix << "\n";
//    std::ifstream dfans_in(dfans_prefix);
//    while (true) {
//        int look = n;
//
//        std::cout << dfans_prefix << " " << n << std::endl;
//        system("pause");
//
//        bool success = augm::read_default_answer(dfans_in, look);
//        std::cout << std::boolalpha << success << "\n";
//
//    }
}
