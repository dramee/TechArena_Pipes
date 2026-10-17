// sim.cpp: локальный симулятор политик укладки труб (без интерактора).
// Сборка: clang++ -O2 -std=c++17 sim.cpp -o sim     (или g++)
// Запуск: ./sim tests/1.in tests/2.in ...           (или ./sim tests/*.in)
// Отладка: clang++ -g -O1 -std=c++17 -fsanitize=address,undefined sim.cpp -o sim_dbg
#include <algorithm>
#include <cstdio>
#include <iterator>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

enum Mode { BEST, FIRST, WORST, SMALL_BIG, PCT_BIG };

struct Cfg {
    string name;
    Mode mode;
    bool right;          // прижимать трубу к правому краю дырки
    double eps;          // SMALL_BIG: "маленькая", если L < eps * stat
    long long switchAt;  // SMALL_BIG: до этого числа труб среднее, дальше медиана
    double pct;          // PCT_BIG: доля прежних труб короче L, с которой берём самую большую дырку
};

static const long long MAX_TRENCH = 100000000LL;
static const int MAXL = 40;

struct Allocator {
    Cfg cfg;
    unordered_map<int, int> byStart, byEnd;  // start -> end, end -> start
    set<pair<int, int>> bySize;              // (размер, start)
    set<int> buckets[MAXL + 2];              // размеры 1..40, корзина 41 = больше 40
    long long lenCount[MAXL + 1];
    long long cnt = 0, sumLen = 0, trench = 0;

    explicit Allocator(const Cfg& c) : cfg(c) {
        for (auto& x : lenCount) x = 0;
        byStart.reserve(1 << 16);
        byEnd.reserve(1 << 16);
    }

    static int bucketOf(int size) { return min(size, MAXL + 1); }

    void addFree(int s, int e) {
        byStart[s] = e;
        byEnd[e] = s;
        buckets[bucketOf(e - s)].insert(s);
        bySize.insert({e - s, s});
    }

    int removeFree(int s) {
        auto it = byStart.find(s);
        if (it == byStart.end()) return s;  // не должно случаться
        int e = it->second;
        byStart.erase(it);
        byEnd.erase(e);
        buckets[bucketOf(e - s)].erase(s);
        bySize.erase(make_pair(e - s, s));
        return e;
    }

    int kth(long long k) const {
        long long acc = 0;
        for (int v = 1; v <= MAXL; v++) {
            acc += lenCount[v];
            if (k < acc) return v;
        }
        return MAXL;
    }
    double median() const { return (kth((cnt - 1) / 2) + kth(cnt / 2)) / 2.0; }
    double mean() const { return (double)sumLen / (double)cnt; }
    double stat() const { return cnt <= cfg.switchAt ? mean() : median(); }

    int worstFit(int L) const {
        if (bySize.empty()) return -1;
        auto it = prev(bySize.end());
        return it->first >= L ? it->second : -1;
    }
    int bestFit(int L) const {
        auto it = bySize.lower_bound(make_pair(L, -1));
        return it == bySize.end() ? -1 : it->second;
    }
    int firstFit(int L) const {
        int s = -1;
        for (int b = L; b <= MAXL + 1; b++)
            if (!buckets[b].empty() && (s < 0 || *buckets[b].begin() < s)) s = *buckets[b].begin();
        return s;
    }

    int pickHole(int L) const {
        switch (cfg.mode) {
            case BEST: return bestFit(L);
            case FIRST: return firstFit(L);
            case WORST: return worstFit(L);
            case SMALL_BIG: {
                if (cnt > 0 && (double)L < cfg.eps * stat()) {
                    int s = worstFit(L);
                    if (s >= 0) return s;
                }
                return bestFit(L);
            }
            case PCT_BIG: {
                if (cnt > 0) {
                    long long below = 0;
                    for (int v = 1; v < L; v++) below += lenCount[v];
                    if ((double)below / (double)cnt >= cfg.pct) {
                        int s = worstFit(L);
                        if (s >= 0) return s;
                    }
                }
                return bestFit(L);
            }
        }
        return -1;
    }

    // возвращает start трубы или -1, если траншея превысила лимит
    long long place(int L) {
        int s = pickHole(L);
        long long D;
        if (s >= 0) {
            int e = removeFree(s);
            if (cfg.right) {
                D = e - L;
                if (s < D) addFree(s, (int)D);
            } else {
                D = s;
                if (s + L < e) addFree(s + L, e);
            }
        } else {
            long long ts = trench;
            auto it = byEnd.find((int)trench);
            if (it != byEnd.end()) {
                ts = it->second;
                removeFree((int)ts);
            }
            long long need = L - (trench - ts);
            long long n = (need + 19) / 20;
            trench += 20 * n;
            if (trench > MAX_TRENCH) return -1;
            D = ts;
            if (ts + L < trench) addFree((int)(ts + L), (int)trench);
        }
        lenCount[L]++;
        sumLen += L;
        cnt++;
        return D;
    }

    void release(int s, int L) {
        int ns = s, ne = s + L;
        auto li = byEnd.find(s);
        if (li != byEnd.end()) {
            int ls = li->second;
            removeFree(ls);
            ns = ls;
        }
        if (byStart.count(ne)) ne = removeFree(ne);
        addFree(ns, ne);
    }
};

// KPI или отрицательное значение при ошибке (причина пишется в why)
static double runPolicy(const vector<int>& cmds, const Cfg& cfg, string& why) {
    Allocator A(cfg);
    vector<int> ps(1, 0), pl(1, 0);
    vector<char> alive(1, 0);
    long long aliveSum = 0;

    for (size_t i = 0; i < cmds.size(); i++) {
        int c = cmds[i];
        if (c == 0) break;
        if (c > 0) {
            if (c > MAXL) { why = "длина трубы > 40 в команде " + to_string(i); return -1; }
            long long D = A.place(c);
            if (D < 0) { why = "траншея длиннее 1e8"; return -1; }
            ps.push_back((int)D);
            pl.push_back(c);
            alive.push_back(1);
            aliveSum += c;
        } else {
            long long k = -(long long)c;
            if (k >= (long long)ps.size()) { why = "удаление PID " + to_string(k) + " ещё не выданной трубы (команда " + to_string(i) + ")"; return -1; }
            if (!alive[k]) { why = "повторное удаление PID " + to_string(k) + " (команда " + to_string(i) + ")"; return -1; }
            alive[k] = 0;
            aliveSum -= pl[k];
            A.release(ps[k], pl[k]);
        }
    }
    if (A.trench == 0) return 0;
    return (double)aliveSum / (double)A.trench;
}

int main(int argc, char** argv) {
    vector<Cfg> cfgs = {
        {"best-left",  BEST,      false, 0,   0,    0},
        {"first-left", FIRST,     false, 0,   0,    0},
        {"worst-left", WORST,     false, 0,   0,    0},
        {"best-right", BEST,      true,  0,   0,    0},
        {"small->big", SMALL_BIG, false, 0.5, 1000, 0},
        {"pct.7->big", PCT_BIG,   false, 0,   0,    0.7},
        {"pct.9->big", PCT_BIG,   false, 0,   0,    0.9},
    };

    if (argc < 2) {
        fprintf(stderr, "usage: %s tests/*.in\n", argv[0]);
        return 1;
    }

    printf("%-22s", "test");
    for (auto& c : cfgs) printf("%-13s", c.name.c_str());
    printf("\n");
    fflush(stdout);

    vector<double> total(cfgs.size(), 0.0);   // сумма KPI
    vector<double> wtotal(cfgs.size(), 0.0);  // сумма KPI * вес
    double wsum = 0.0;                        // сумма весов
    int files = 0;

    for (int i = 1; i < argc; i++) {
        FILE* f = fopen(argv[i], "r");
        if (!f) { fprintf(stderr, "не удалось открыть %s\n", argv[i]); continue; }
        vector<int> cmds;
        int x;

        // Первая строка открытых тестов: "<a> <b>" (например "1 9997"), это заголовок, а не команды.
        // Если в первой строке два числа, пропускаем её; если одно, считаем его командой.
        // Формат заголовка: "<вес теста, может быть дробным> <число команд>", например "0.1 6" или "1 9997".
        char line[256];
        double weight = 1.0;
        long long hdrB = -1;
        if (fgets(line, sizeof line, f)) {
            double a;
            long long b;
            int got = sscanf(line, "%lf %lld", &a, &b);
            if (got == 2) { weight = a; hdrB = b; }
            else if (got == 1) cmds.push_back((int)a);   // нет заголовка: одно число = команда
        }
        while (fscanf(f, "%d", &x) == 1) {
            cmds.push_back(x);
            if (x == 0) break;
        }
        fclose(f);
        // hdrB, скорее всего, число команд: предупреждаем, если не сходится (это только предупреждение)
        if (hdrB >= 0 && (long long)cmds.size() != hdrB && (long long)cmds.size() != hdrB + 1)
            fprintf(stderr, "\n[%s] заголовок: вес %g, команд %lld, прочитано: %zu\n", argv[i], weight, hdrB, cmds.size());

        string name = argv[i];
        size_t p = name.find_last_of("/\\");
        if (p != string::npos) name = name.substr(p + 1);
        printf("%-22s", name.c_str());
        fflush(stdout);  // видно, на каком тесте остановилось

        for (size_t j = 0; j < cfgs.size(); j++) {
            string why;
            double kpi = runPolicy(cmds, cfgs[j], why);
            if (kpi < 0) {
                printf("%-13s", "INVALID");
                fprintf(stderr, "\n[%s / %s] %s\n", name.c_str(), cfgs[j].name.c_str(), why.c_str());
            } else {
                printf("%-13.4f", kpi);
                total[j] += kpi;
                wtotal[j] += kpi * weight;   // INVALID даёт KPI 0
            }
            fflush(stdout);
        }
        printf("  w=%g\n", weight);
        wsum += weight;
        files++;
    }

    printf("%-22s", "mean");
    for (size_t j = 0; j < cfgs.size(); j++) printf("%-13.4f", files ? total[j] / files : 0.0);
    printf("\n");
    printf("%-22s", "weighted mean");
    for (size_t j = 0; j < cfgs.size(); j++) printf("%-13.4f", wsum > 0 ? wtotal[j] / wsum : 0.0);
    printf("\n");
    printf("%-22s", "score (sum KPI*w)");
    for (size_t j = 0; j < cfgs.size(); j++) printf("%-13.4f", wtotal[j]);
    printf("\n");
    return 0;
}