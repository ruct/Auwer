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
    void write() {
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

const int GENS = 200;
const int REAP = 4;
const int PEACE = 1000;
chromo scent[PEACE];
const int ST_FEW_SWAPS = 50;
const int FSKIP = 5;
const int TSKIP = 5;
namespace GA {

    void init() {
        for (int i = 0; i < PEACE; ++i) {
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

//    int del_rnd(std::vector <int>& v) {
//        int i = gen()%v.size();
//        std::swap(v[i], v.back());
//        i = v.back();
//        v.pop_back();
//        return i;
//    }

    bool WRITE_CROSS = false;
    int CUR_GEN_CNT_FIRST_PHASE;
    int CUR_GEN_CNT_SECOND_PHASE;
    int CUR_GEN_CNT_FIRST_MUTATION;
    int CUR_GEN_CNT_THIRD_MUTATION;
    int OVERALL_GEN;
    int CUR_GEN;

    void cross(chromo& a, chromo& b, chromo& t) {
        temp_chromo temp;


        if (WRITE_CROSS) {
            std::cout << "CROSS: \n";
            std::cout << "A: \n"; a.write();
            std::cout << "B: \n"; b.write();
        }

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
                ++OVERALL_GEN;
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
//                lft.erase(std::find(lft.begin(), lft.end(), id));
//                std::cout << std::boolalpha << happened << " LDD" << std::endl;
//                std::cout << "T" << " " << i << " " << j << std::endl;
//                std::cout << bound.size() << std::endl;
                bound.emplace_back(i, j);
//                std::cout << bound[0].first << " " << bound[0].second << std::endl;
//                std::cout << "G" << std::endl;

//                std::cout << bound[0].first << " " << bound[0].second << std::endl;
                bound_update();
//                std::cout << "new size = " << bound.size() << std::endl;

//                std::cout << "M" << std::endl;
            };

        place_fragment(cntsq, cntsq, gen()%(cntsq*cntsq));

        while (placed < cntsq*cntsq) {
//            if (WRITE_CROSS && placed % 3 == 0)
//                system("cls");

            if (WRITE_CROSS) {
                std::cout << "AA:\n";
                a.write();
                std::cout << "\n";
                std::cout << "BB:\n";
                b.write();
                std::cout << "\n";
                temp.write();
                system("pause");
                std::cout << "Placed " << placed << std::endl;
            }

            //{index in number, direction, suitable_id}
            std::vector <std::tuple <int, int, int> > suitable;
            ///  first phase
            if (CUR_GEN > 1){
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
                            if (WRITE_CROSS) {
                                std::cout << "Excellent match " << cur_id << " " << a.perm[ni1][nj1] << " " << d << std::endl;
                            }
                        }
                    }
                }

                if (WRITE_CROSS)
                    std::cout << "First phase " << suitable.size() << std::endl;

                if (suitable.size()) {
                    if (gen()%100 <= FSKIP) {
                        int dirs[] = {0, 1, 2, 3};
                        std::shuffle(dirs, dirs+4, gen);

                        int P = gen()%bound.size();
                        int i = bound[P].first,
                            j = bound[P].second;
                        for (int _d = 0; _d < 4; ++_d) {
                            int d = dirs[_d];
                            int ni = i+dx[d],
                                nj = j+dy[d];
                            if (!temp.can(ni, nj))
                                continue;
                            ++CUR_GEN_CNT_FIRST_MUTATION;
                            place_fragment(ni, nj, lft[gen()%lft.size()]);
                        }
                        continue;
                    }

                    int srnd = gen()%suitable.size();
                    int TT = std::get<0> (suitable[srnd]);
                    int id = std::get<2> (suitable[srnd]);
                    int d = std::get<1> (suitable[srnd]);

                    int i = bound[TT].first,
                        j = bound[TT].second;

                    if (WRITE_CROSS) {
                        system("pause");
                    }

                    ++CUR_GEN_CNT_FIRST_PHASE;
                    place_fragment(i+dx[d], j+dy[d], id);

                    continue;
                }
            }

            ///second phase
            {
                suitable.clear();
                //  best {direction, TT}
                std::vector <int> best_for_bb(4*cntsq*cntsq);

//                std::cout << "START" << std::endl;
                for (int BB : lft) {
                    ld mn_for_dirs[4] = {1e9, 1e9, 1e9, 1e9};
                    int mn_for_dirs_id[4] = {-1, -1, -1, -1};

                    for (int TT = 0; TT < bound.size(); ++TT) {
                        int i = bound[TT].first,
                            j = bound[TT].second;

                        for (int d = 0; d < 4; ++d) {
                            if (!temp.can(i+dx[d], j+dy[d]))
                                continue;
                            ld cur_match = pair_match(BB, temp.perm[i][j], d^1);
                            if (cur_match < mn_for_dirs[d]) {
                                mn_for_dirs[d] = cur_match;
                                mn_for_dirs_id[d] = TT;
                            }
                        }
                    }
                    for (int d = 0; d < 4; ++d)
                        best_for_bb[4*BB+d] = mn_for_dirs_id[d];
                }
//                std::cout << "MIDDLE" << std::endl;
                for (int TT = 0; TT < bound.size(); ++TT) {
                    int i = bound[TT].first, j = bound[TT].second;

                    int cur_id = temp.perm[i][j];
                    ld mn_for_dirs[4] = {1e9, 1e9, 1e9, 1e9};
                    int mn_for_dirs_id[4] = {-1, -1, -1, -1};

//                    std::cout << "finding best BB" << std::endl;
                    for (int BB = 0; BB < lft.size(); ++BB) {
                        for (int d = 0; d < 4; ++d) {
                            if (!temp.can(i+dx[d], j+dy[d]))
                                continue;
                            ld cur_match = pair_match(cur_id, lft[BB], d);
                            if (cur_match < mn_for_dirs[d]) {
                                mn_for_dirs[d] = cur_match;
                                mn_for_dirs_id[d] = BB;
                            }
                        }
                    }
//                    std::cout << "found best BB " << std::endl;


                    int pos1 = a.wher[cur_id],
                        pos2 = b.wher[cur_id];

                    for (int d = 0; d < 4; ++d) {
                        if (!temp.can(i+dx[d], j+dy[d]))
                            continue;

                        int ni1 = (pos1/cntsq)+dx[d], nj1 = (pos1%cntsq)+dy[d];
                        int ni2 = (pos2/cntsq)+dx[d], nj2 = (pos2%cntsq)+dy[d];

                        if (std::min(ni1, nj1) >= 0 && std::max(ni1, nj1) < cntsq) {
                            int id1 = a.perm[ni1][nj1];
                            if (!used[id1] && mn_for_dirs_id[d] == id1 && best_for_bb[4*id1+d] == TT) {
                                suitable.emplace_back(TT, d, id1);
                            }
                        }
                        if (std::min(ni2, nj2) >= 0 && std::max(ni2, nj2) < cntsq) {
                            int id2 = b.perm[ni2][nj2];
                            if (!used[id2] && mn_for_dirs_id[d] == id2 && best_for_bb[4*id2+d] == TT) {
                                suitable.emplace_back(TT, d, id2);
                            }
                        }
                    }
                }
//                std::cout << "FINISH" << std::endl;

//                std::cout << "Second phase " << suitable.size() << std::endl;
                if (suitable.size()) {
                    int srnd = gen()%suitable.size();
                    int TT = std::get<0> (suitable[srnd]);
                    int id = std::get<2> (suitable[srnd]);
                    int d = std::get<1> (suitable[srnd]);

                    int i = bound[TT].first,
                        j = bound[TT].second;

                    if (WRITE_CROSS) {
                        system("pause");
                    }

                    ++CUR_GEN_CNT_SECOND_PHASE;
                    place_fragment(i+dx[d], j+dy[d], id);
                    continue;
                }
            }

            ///third phase
            {
                if (WRITE_CROSS) {
                    std::cout << "third - start" << std::endl;
                    std::cout << bound.size() << std::endl;
                }

                if (gen()%100 <= TSKIP) {
                    int dirs[] = {0, 1, 2, 3};
                    std::shuffle(dirs, dirs+4, gen);

                    int P = gen()%bound.size();
                    int i = bound[P].first,
                        j = bound[P].second;
                    for (int _d = 0; _d < 4; ++_d) {
                        int d = dirs[_d];
                        int ni = i+dx[d],
                            nj = j+dy[d];
                        if (!temp.can(ni, nj))
                            continue;
                        ++CUR_GEN_CNT_THIRD_MUTATION;
                        place_fragment(ni, nj, lft[gen()%lft.size()]);
                    }
                    continue;
                }

                int TT = gen()%bound.size();
                int i = bound[TT].first,
                    j = bound[TT].second;
                int id = temp.perm[i][j];

                if (WRITE_CROSS)
                    std::cout << "Third phase " << i << " " << j << std::endl;

                ld mn = 1e9;
                int mn_id = -1, mn_d = -1;
                for (int d = 0; d < 4; ++d) {
                    int ni = i+dx[d],
                        nj = j+dy[d];

                    if (!temp.can(ni, nj))
                        continue;

                    if (WRITE_CROSS)
                        std::cout << "random_find " << id << " " << d << " " << ni << " " << nj << std::endl;

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
        t.eval = -1;
        t.mk_wher();

        ld F1 = fit_cost(a), F2 = fit_cost(b);
//        std::cout << "CROSSSS " << F1 << " " << F2 << " " << fit_cost(t) << std::endl;
        a.eval = b.eval = -1;
        ld S1 = fit_cost(a), S2 = fit_cost(b);

        if (F1 != S1 || F2 != S2) {
            exit(-1);
        }
    }

    void breed() {
//        std::cout << "start breeding" << std::endl;

        int size = 0;
        chromo cur[REAP];
        for (int i = 0; i < REAP; ++i)
            cur[i] = scent[i];

        while (size < PEACE) {
            for (int i = 0; size < PEACE && i < REAP; ++i)
            for (int j = 0; size < PEACE && j < REAP; ++j) {
                if (i == j) continue;
                cross(cur[i], cur[j], scent[size]);
                ++size;
            }
        }

//        while (size < (PEACE)/2) {
//            for (int i = 0; i < REAP; i += 2) {
//                cross(scent[i], scent[i+1], scent[size]);
//                ++size;
//                if (size == PEACE)
//                    break;
//            }
////            std::cout << size << std::endl;
//        }
        while (size < PEACE) {
            int i = gen()%size,
                j = gen()%(size-1);
            if (j >= i) ++j;
            cross(scent[i], scent[j], scent[size]);
            ++size;
        }

        std::shuffle(scent, scent+PEACE, gen);
//        std::cout << "breed ended" << std::endl;
//        for (int i = 0; i < PEACE; ++i)
//            scent[i] = scent[i+REAP];
    }
    void harvest() {
//        std::cout << "start harvesting" << std::endl;

        std::vector <std::pair <ld, int> > _scent(PEACE);
        for (int i = 0; i < PEACE; ++i)
            _scent[i] = {fit_cost(scent[i]), i};
        std::sort(_scent.begin(), _scent.end());
        for (int i = 0; i < REAP; ++i)
            std::cout << _scent[i].first << " ";
        std::cout << std::endl;
//        scent[0].write();

        std::vector <chromo> buf(REAP);
        for (int i = 0; i < REAP; ++i)
            buf[i] = scent[_scent[i].second];
        for (int i = 0; i < REAP; ++i)
            scent[i] = buf[i];

//        std::cout << "harvested" << std::endl;
    }

    int16_t GA_ans[cntsq][cntsq];
    void genetic_algorithm() {
        init();
        harvest();

        std::cout << fit_cost(scent[0]) << std::endl;
        system("pause");

//        WRITE_CROSS = 1;
        std::cout << "we are here" << std::endl;
        for (CUR_GEN = 0; CUR_GEN < GENS; ++CUR_GEN) {
            breed(), harvest();
            std::cout << CUR_GEN_CNT_FIRST_PHASE << " " << OVERALL_GEN << " "
                << double(CUR_GEN_CNT_FIRST_PHASE)/double(OVERALL_GEN) << std::endl;
            std::cout << CUR_GEN_CNT_SECOND_PHASE << " " << OVERALL_GEN << " "
                << double(CUR_GEN_CNT_SECOND_PHASE)/double(OVERALL_GEN) << std::endl;
            std::cout << CUR_GEN_CNT_FIRST_MUTATION << " " << OVERALL_GEN << " "
                << double(CUR_GEN_CNT_FIRST_MUTATION)/double(OVERALL_GEN) << std::endl;
            std::cout << CUR_GEN_CNT_THIRD_MUTATION << " " << OVERALL_GEN << " "
                << double(CUR_GEN_CNT_THIRD_MUTATION)/double(OVERALL_GEN) << std::endl;

            CUR_GEN_CNT_FIRST_PHASE = 0;
            CUR_GEN_CNT_SECOND_PHASE = 0;
            CUR_GEN_CNT_FIRST_MUTATION = 0;
            CUR_GEN_CNT_THIRD_MUTATION = 0;
            OVERALL_GEN = 0;
//            if (i == 25) {
//                system("cls");
//                system("pause");
//                system("pause");
//                scent[0].write();
//                std::cout << "\n";
//                scent[1].write();
//
//                for (int i = 0; i < cntsq*cntsq; ++i)
//                    std::cout << scent[0].wher[i] << " ";
//                std::cout << "\n";
//                for (int i = 0; i < cntsq*cntsq; ++i)
//                    std::cout << scent[1].wher[i] << " ";
//                std::cout << "\n";
//
//                std::cout << "\nChildren:\nfirst:\n";
//                chromo cur;
//
//                WRITE_CROSS = 1;
//                cross(scent[0], scent[1], cur);
//                cur.write();
//                system("pause");
//                system("pause");
//                system("pause");
//
//                std::cout << "second:\n";
//                cross(scent[0], scent[1], cur);
//                cur.write();
//                WRITE_CROSS = 0;
//                std::cout << "END" << std::endl;
//                system("pause");
//            }
            ld sum = 0, mn = 1e18;
            for (int j = 0; j < REAP; ++j) {
                sum += fit_cost(scent[j]);
                mn = std::min(mn, fit_cost(scent[j]));
            }
            std::cout << "Generation #" << CUR_GEN << " " << sum << " " << mn << std::endl;
//            system("pause");

        }

        std::cout << "\n\n\n";
        scent[0].write();
        for (int i = 0; i < cntsq; ++i)
        for (int j = 0; j < cntsq; ++j)
            GA_ans[i][j] = scent[0].perm[i][j];
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

    int n = 1200;
    string shimg = shimg_prefix+"\\"+std::to_string(n)+".txt";
    string crimg = crimg_prefix+"\\"+std::to_string(n)+".txt";
    std::ifstream shimg_in(shimg, std::ios::binary);
    std::ifstream crimg_in(crimg);

    augm::read_picture(shimg_in);
    GA::recalc();

    GA::genetic_algorithm();
    string ans_write_path = R"(C:\Users\Main\Base\huawei\parsed\)" + NAME + "\\" +
            std::to_string(size) + "\\answers.txt";
    std::ofstream ans_out(ans_write_path);
    ans_out << std::to_string(n) << ".png" << std::endl;
    for (int i = 0; i < cntsq; ++i)
    for (int j = 0; j < cntsq; ++j)
        ans_out << GA::GA_ans[i][j] << " ";
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
