#include <tuple>
#include <future>
#include <chrono>
#include <random>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>

using std::vector;
using std::string;
using std::cout;
using std::endl;

std::mt19937 GGEN(3711+time(NULL));
const int WH = 512;
const int SIZE = 16;
const string NAME = "data_test1_blank";
const int NL = 1810, NR = 2099;
const int CNTSQ = WH/SIZE;
typedef double ld;

namespace augm {
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

    int picture_buf[3][WH][WH];
    void read_picture(std::ifstream& in) {
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
            result += ld(sqr(abs(a[c][T]-b[c][T])));
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
        for (int d = 0; d < 4; ++d)
        for (int p = 0; p < CNTSQ*CNTSQ; ++p) {
            ld *b = diss[d][p] + (p == 0);
            for (int j = 0; j < CNTSQ*CNTSQ; ++j)
                if (p != j && diss[d][p][j] < *b)
                    b = diss[d][p]+j;
            bbuddy[d][p] = b-diss[d][p];
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
        void write() const {
            for (int i = 0; i < CNTSQ; ++i) {
                for (int j = 0; j < CNTSQ; ++j)
                    cout << perm[i][j] << " ";
                cout << endl;
            }
        }
        void write(std::ofstream& out, string s1 = "\t", string s2 = "\n") const {
            for (int i = 0; i < CNTSQ; ++i) {
                for (int j = 0; j < CNTSQ; ++j)
                    out << perm[i][j] << s1;
                out << s2;
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
            return eval;
        }
    };

    void random_chromo(chromo& trg) {
        vector <int> cont(CNTSQ*CNTSQ);
        std::iota(cont.begin(), cont.end(), 0);
        std::shuffle(cont.begin(), cont.end(), GGEN);

        for (int i = 0; i < CNTSQ; ++i)
        for (int j = 0; j < CNTSQ; ++j)
            trg.perm[i][j] = cont[i*CNTSQ+j];
        trg.eval = -1;
        trg.make_wher();
    }
}


const int CORES = 8;
const int GENS = 40;
const int WAVES = 1;
const int POP = 150;
const int FSINCE = 0;
const int PREAP = 20;
const int MRATE = 5;
const int ODDS = 5;

//  population & kernel
namespace ga {
    struct population {
        int nogen;
        chromo emin;
        vector <chromo> stock;
        vector <chromo> spring;
        void init(bool emin_reset = true) {
            nogen = 0;
            stock.resize(POP);

            for (int i = 0; i < POP; ++i)
                random_chromo(stock[i]);
            if (emin_reset)
                emin = stock[0];
        }
    };

    const int TCNTSQ = 2*CNTSQ+1;
    struct kernel {

    private:
        vector <int> edge, plan;
        int perm[TCNTSQ][TCNTSQ];
        int mni = TCNTSQ, mxi = -1;
        int mnj = TCNTSQ, mxj = -1;
        bool was[CNTSQ*CNTSQ];

    public:
        kernel();
        void write();
        void convert(chromo&);

        bool used(int id) const { return was[id]; }
        vector <int> const& get_edge() const { return edge; }
        vector <int> const& get_plan() const { return plan; }
        int arranged() const { return CNTSQ*CNTSQ-int(plan.size()); }

        //  assuming i in [0, TCNTSQ), j in [0, TCNTSQ)
        int get(int i, int j) const { return perm[i][j]; }
        int can(int i, int j) const {
            if (perm[i][j] != -1)
                return false;

            if (mxi-mni+1 == CNTSQ && (i < mni || i > mxi))
                return false;
            if (mxj-mnj+1 == CNTSQ && (j < mnj || j > mxj))
                return false;
            return true;
        }
        void place(int i, int j, int t) {
            assert(perm[i][j] == -1);

            perm[i][j] = t;
            mni = std::min(mni, i);
            mnj = std::min(mnj, j);
            mxi = std::max(mxi, i);
            mxj = std::max(mxj, j);
            edge.push_back(i*TCNTSQ+j);

            plan.erase(std::find(plan.begin(), plan.end(), t));

            vector <int> new_edge;
            for (int pos : edge) {
                int ei = pos/TCNTSQ,
                    ej = pos%TCNTSQ;
                int cnt = can(ei+1, ej)+can(ei-1, ej)
                        + can(ei, ej+1)+can(ei, ej-1);

                if (cnt != 0)
                    new_edge.push_back(ei*TCNTSQ+ej);
            }
            edge.swap(new_edge);
            was[t] = true;
        }
    };
    kernel::kernel() {
        for (int i = 0; i < TCNTSQ; ++i)
        for (int j = 0; j < TCNTSQ; ++j)
            perm[i][j] = -1;
        plan.resize(CNTSQ*CNTSQ);
        std::iota(plan.begin(), plan.end(), 0);
    }
    void kernel::write() {
        for (int i = 0; i < TCNTSQ; ++i) {
            for (int j = 0; j < TCNTSQ; ++j)
                cout << perm[i][j] << "\t";
            cout << endl;
        }
    }
    void kernel::convert(chromo& t) {
        for (int i = 0; i < CNTSQ; ++i)
        for (int j = 0; j < CNTSQ; ++j)
            t.perm[i][j] = perm[i+mni][j+mnj];
        t.make_wher();
        t.eval = -1;
    }

}

//  chromosome shifting
namespace ga {
    void shiftR(chromo& t, int i) {
        int buf = t.perm[i][CNTSQ-1];
        for (int j = CNTSQ-1; j != 0; --j)
            t.perm[i][j] = t.perm[i][j-1];
        t.perm[i][0] = buf;
    }
    void shiftD(chromo& t, int j) {
        int buf = t.perm[CNTSQ-1][j];
        for (int i = CNTSQ-1; i != 0; --i)
            t.perm[i][j] = t.perm[i-1][j];
        t.perm[0][j] = buf;
    }

    void shiftL(chromo& t, int i) {
        int buf = t.perm[i][0];
        for (int j = 0; j+1 < CNTSQ; ++j)
            t.perm[i][j] = t.perm[i][j+1];
        t.perm[i][CNTSQ-1] = buf;
    }
    void shiftU(chromo& t, int j) {
        int buf = t.perm[0][j];
        for (int i = 0; i+1 < CNTSQ; ++i)
            t.perm[i][j] = t.perm[i+1][j];
        t.perm[CNTSQ-1][j] = buf;
    }
    void shift(chromo& t, int d) {
        int i = GGEN()%CNTSQ;
        switch (d) {
            case 0 : shiftR(t, i);
            case 1 : shiftL(t, i);
            case 2 : shiftD(t, i);
            case 3 : shiftU(t, i);
        }
        t.make_wher();
        t.eval = -1;
    }
}


ld fphase_time;
ld sphase_time;
ld tphase_time;

//  genetic operators
namespace ga {
    void print(std::ofstream& out, int nofimage, chromo const& t) {
        out << nofimage << ".png\n";
        t.write(out, " ", "");
        out << endl;
    }

    void mutate(chromo& t, std::mt19937& gen) {
    //  Vertical shift chance
    static const int VSHIFT = 50;

        if (gen()%100 < VSHIFT)
            shift(t, gen()&1);
        else
            shift(t, 2^(gen()&1));
    }

    bool FLAG = 0;
    bool first_phase(population const& pop, chromo const& a, chromo const& b, kernel& step, std::mt19937& gen) {
        augm::TimeLogger fphase_logger(&fphase_time);

        if (FLAG) {
            a.write(), cout << "\n", b.write();
            cout << endl << endl;
            step.write();
            cout << endl;
            system("pause");
        }

        if (pop.nogen <= FSINCE)
            return false;
        auto& edge = step.get_edge();
        auto& plan = step.get_plan();

        //  {target_place(i, j), target_id}
        vector <std::tuple <int, int, int> > fsuit;
        for (int T = 0; T < edge.size(); ++T) {
            int pT = edge[T];
            int i = pT/TCNTSQ, j = pT%TCNTSQ;

            int id = step.get(i, j);
            int p1 = a.wher[id], p2 = b.wher[id];
            int i1 = p1/CNTSQ, j1 = p1%CNTSQ;
            int i2 = p2/CNTSQ, j2 = p2%CNTSQ;
            for (int d = 0; d < 4; ++d) {
                if (!step.can(i+dx[d], j+dy[d]))
                    continue;

                int ni1 = i1+dx[d], nj1 = j1+dy[d];
                if (std::min(ni1, nj1) < 0 ||
                    std::max(ni1, nj1) >= CNTSQ)
                    continue;

                int ni2 = i2+dx[d], nj2 = j2+dy[d];
                if (std::min(ni2, nj2) < 0 ||
                    std::max(ni2, nj2) >= CNTSQ)
                    continue;

                if (!step.used(a.perm[ni1][nj1]) && a.perm[ni1][nj1] == b.perm[ni2][nj2])
                    fsuit.emplace_back(i+dx[d], j+dy[d], a.perm[ni1][nj1]);
            }
        }

        if (fsuit.empty())
            return false;

        if (FLAG) {
            cout << "writing fsuit: \n";
            for (auto i : fsuit)
                cout << std::get<0> (i) << " " << std::get<1> (i) << " " << std::get<2> (i) << endl;
            system("pause");
        }

        int pi, pj, pid;
        std::tie(pi, pj, pid) = fsuit[gen()%fsuit.size()];
        return step.place(pi, pj, pid), true;
    }
    bool second_phase(population const& pop, chromo const& a, chromo const& b, kernel& step, std::mt19937& gen) {
        augm::TimeLogger sphase_logger(&sphase_time);

        if (FLAG)
            cout << "SECOND PHASE" << endl;

        auto& edge = step.get_edge();

        //  {target_place(i, j), target_id}
        vector <std::tuple <int, int, int> > ssuit;
        for (int T = 0; T < edge.size(); ++T) {
            int pT = edge[T];
            int i = pT/TCNTSQ, j = pT%TCNTSQ;

            int id = step.get(i, j);
            int p1 = a.wher[id], p2 = b.wher[id];
            int i1 = p1/CNTSQ, j1 = p1%CNTSQ;
            int i2 = p2/CNTSQ, j2 = p2%CNTSQ;
            for (int d = 0; d < 4; ++d) {
                if (!step.can(i+dx[d], j+dy[d]))
                    continue;

                int ni1 = i1+dx[d], nj1 = j1+dy[d];
                if (std::min(ni1, nj1) >= 0 && std::max(ni1, nj1) < CNTSQ) {
                    int pid = a.perm[ni1][nj1];

                    if (!step.used(pid) && bbuddy[d][id] == pid && bbuddy[d^1][pid] == id)
                        ssuit.emplace_back(i+dx[d], j+dy[d], pid);
                }

                int ni2 = i2+dx[d], nj2 = j2+dy[d];
                if (std::min(ni2, nj2) >= 0 && std::max(ni2, nj2) < CNTSQ) {
                    int pid = b.perm[ni2][nj2];
                    if (!step.used(pid) && bbuddy[d][id] == pid && bbuddy[d^1][pid] == id)
                        ssuit.emplace_back(i+dx[d], j+dy[d], pid);
                }
            }
        }

        if (ssuit.empty())
            return false;

        int pi, pj, pid;
        std::tie(pi, pj, pid) = ssuit[gen()%ssuit.size()];
        return step.place(pi, pj, pid), true;
    }
    bool third_phase(population const& pop, chromo const& a, chromo const& b, kernel& step, std::mt19937& gen) {
        augm::TimeLogger tphase_logger(&tphase_time);

        if (FLAG)
            cout << "THIRD PHASE" << endl;

        auto &edge = step.get_edge(),
            &plan = step.get_plan();

        int pT = edge[gen()%edge.size()];
        int i = pT/TCNTSQ, j = pT%TCNTSQ;
        int id = step.get(i, j);

        std::tuple <ld, int, int, int> tsuit(1e15, -1, -1, -1);
        for (int d = 0; d < 4; ++d) {
            int ni = i+dx[d], nj = j+dy[d];
            if (!step.can(ni, nj))
                continue;
            for (int pid : plan) {
                ld cur_diss = diss[d][id][pid];
                auto ctuple = std::make_tuple(cur_diss, ni, nj, pid);
                tsuit = std::min(tsuit, ctuple);
            }
        }

        ld mn_diss;
        int pi, pj, pid;
        std::tie(mn_diss, pi, pj, pid) = tsuit;
        return step.place(pi, pj, pid), true;
    }

    int FPHASE = 0;
    int SPHASE = 0;
    int TPHASE = 0;
    int APHASE = 0;
    std::ofstream* for_check;
    void cross(population const& pop, chromo const& a, chromo const& b, chromo& t, std::mt19937& gen) {
        kernel step;
        step.place(CNTSQ, CNTSQ, gen()%(CNTSQ*CNTSQ));

        while (step.arranged() != CNTSQ*CNTSQ) {
            if (FLAG) {
                a.write(), cout << "\n", b.write();
                step.write();
                cout << "\n\n";
                system("pause");
            }

            ++APHASE;
            if (first_phase(pop, a, b, step, gen)) {
                ++FPHASE;
                continue;
            }
            if (second_phase(pop, a, b, step, gen)) {
                ++SPHASE;
                continue;
            }
            if (third_phase(pop, a, b, step, gen)) {
                ++TPHASE;
                continue;
            }

            assert(false);
        }
        step.convert(t);
    }

    void preed(population& pop, int l, int r) {
        std::mt19937 gen; {
            auto t = std::chrono::high_resolution_clock::now();
            gen = std::mt19937(time(NULL)+t.time_since_epoch().count());
        }
        vector <chromo> const& stock = pop.stock;
        for (int i = l; i <= r; ++i) {
            int a = gen()%stock.size(),
                b = gen()%(int(stock.size())-1);
            if (b >= a) ++b;
            cross(pop, stock[a], stock[b], pop.spring[i], gen);
        }
    }
    void breed(population& pop) {
        ++pop.nogen;

        FPHASE = 0;
        SPHASE = 0;
        TPHASE = 0;
        APHASE = 0;

        pop.spring.resize(POP-pop.stock.size());

        std::vector <int> ends;
        int each = (POP-pop.stock.size())/CORES;
        for (int i = each; i < POP-pop.stock.size(); i += each)
            ends.push_back(i);
        ends.back() = POP-pop.stock.size();

        std::future <void> forks[CORES];
        for (int i = 0; i < CORES; ++i) {
            int l = (i ? ends[i-1] : 0), r = ends[i]-1;
            forks[i] = std::async(std::launch::async, preed, std::ref(pop), l, r);
        }
        for (int i = 0; i < CORES; ++i)
            forks[i].get();

        pop.spring.swap(pop.stock);
        while (pop.spring.size()) {
            pop.stock.push_back(pop.spring.back());
            pop.spring.pop_back();
        }

        for (chromo const& T : pop.stock)
            if (T.fit() < pop.emin.fit())
                pop.emin = T;
    }

    void truncate(population& pop, std::mt19937& gen) {
        ld average = 0;
        for (chromo& t : pop.stock)
            average += t.fit();
        average /= pop.stock.size();
//        cout << "AVERAGE " << average << " " << pop.stock.size() << endl;

        auto& stock = pop.stock;
        std::sort(stock.begin(), stock.end(),
            [](chromo const& a, chromo const& b) {
                return a.fit() < b.fit();
            });

        int odds = (ODDS*POP)/100;
        int retain = (PREAP*POP)/100;

        for (int T = 0; T < odds; ++T) {
            int i = retain+gen()%(POP-retain);
            std::swap(stock[retain], stock[i]);
            ++retain;
        }

        stock.resize(retain);
        for (chromo& t : stock)
            if (gen()%100 < MRATE)
                mutate(t, gen);
    }

    chromo overall(int nofimage = -1) {
        population pop;
//        cout << "initiating" << endl;
        pop.init();

        for (int E = 0; E < WAVES; ++E) {
            for (int T = 0; T < GENS; ++T) {
//                cout << nofimage << " Generation " << T << endl;
                truncate(pop, GGEN);
                breed(pop);
//                cout << "Rate MDA: " << ld(FPHASE)/ld(APHASE) << " " << ld(SPHASE)/ld(APHASE) << " " << ld(TPHASE)/ld(APHASE) << std::endl;

//                cout << "fit_min " << pop.emin.fit() << endl;
            }
            pop.init(false);
        }
        return pop.emin;
    }
}

//  picture & ans paths, etc...
namespace augm {
    const string PATH = R"(C:\Users\Main\Base\huawei\parsed\)";
    const string AUTH = PATH+NAME+"\\"+std::to_string(SIZE);
    const string MINE = PATH+NAME+"\\"+std::to_string(SIZE)+"\\answers.txt";

    std::ofstream mine(MINE);
    void solve(int n) {
        cout << "solving " << n << endl;

        string source = AUTH+"\\"+std::to_string(n)+".txt";
        std::ifstream image(source, std::ios::binary);
        read_picture(image);
        ga::reset(), ga::compute();

        ga::for_check = &mine;
        ga::print(mine, n, ga::overall(n));
        cout << ld(clock())/ld(CLOCKS_PER_SEC) << " " << fphase_time << " " << sphase_time << " " << tphase_time << endl;
    }
}
int main(){
    std::ios::sync_with_stdio(0);
    for (int n = NL; n <= NR; ++n)
        augm::solve(n);
}
