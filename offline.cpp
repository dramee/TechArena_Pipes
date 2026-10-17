// offline.cpp: офлайн-анализ задачи о прокладке трубопровода (без интерактора).
// Читает тест целиком (заголовок "<вес> <число команд>", затем команды), знает будущее и ищет
// раскладку с минимальной итоговой длиной траншеи (KPI = сумма живых труб / траншея).
//
// Что считается:
//   LB      нижняя граница траншеи: пик суммарной длины одновременно лежащих труб, округлённый вверх до 20 м
//   online  обычный best-fit без знания будущего (для сравнения)
//   oracle  best-fit, знающий моменты удаления труб (выбирает дырку и сторону по близости времён удаления соседей)
//   exact   полный перебор для небольших тестов: доказывает оптимум (OPT) или честно говорит, что не доказал
//
// Сборка: clang++ -O2 -std=c++17 offline.cpp -o offline     (или g++)
// Запуск: ./offline [--dump dir] [--w2 X] [--k N] [--max-exact-events N] [--node-limit N] tests/*.in
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;
typedef long long ll;

static const ll MAX_TRENCH = 100000000LL;
static const int MAXL = 40;

struct Pipe { int len, laid, rem; };   // rem: номер события удаления, либо m, если труба остаётся до конца

struct Test {
    string name;
    double weight = 1.0;
    long long hdrCount = -1;
    vector<int> cmds;       // команды без завершающего 0
    int m = 0;
    vector<Pipe> pipes;     // индекс = PID (с 1)
    vector<int> evPid;      // для события-укладки: PID
    ll aliveSum = 0, peak = 0;
    string err;
};

static bool loadTest(const char* path, Test& t) {
    FILE* f = fopen(path, "r");
    if (!f) { t.err = "не удалось открыть файл"; return false; }
    string nm = path;
    size_t p = nm.find_last_of("/\\");
    t.name = p == string::npos ? nm : nm.substr(p + 1);

    char line[256];
    vector<int> raw;
    if (fgets(line, sizeof line, f)) {
        double a; long long b;
        int got = sscanf(line, "%lf %lld", &a, &b);
        if (got == 2) { t.weight = a; t.hdrCount = b; }
        else if (got == 1) raw.push_back((int)a);
    }
    int x;
    while (fscanf(f, "%d", &x) == 1) { raw.push_back(x); if (x == 0) break; }
    fclose(f);
    // обрезаем по первому 0
    for (size_t i = 0; i < raw.size(); i++) { if (raw[i] == 0) break; t.cmds.push_back(raw[i]); }
    t.m = (int)t.cmds.size();

    t.pipes.assign(1, Pipe{0, -1, -1});
    t.evPid.assign(t.m, 0);
    ll live = 0;
    for (int e = 0; e < t.m; e++) {
        int c = t.cmds[e];
        if (c > 0) {
            if (c > MAXL) { t.err = "длина трубы > 40 в команде " + to_string(e); return false; }
            t.pipes.push_back(Pipe{c, e, -1});
            t.evPid[e] = (int)t.pipes.size() - 1;
            live += c;
            t.peak = max(t.peak, live);
        } else {
            ll k = -(ll)c;
            if (k >= (ll)t.pipes.size()) { t.err = "удаление ещё не выданной трубы, команда " + to_string(e); return false; }
            if (t.pipes[k].rem != -1) { t.err = "повторное удаление PID " + to_string(k); return false; }
            t.pipes[k].rem = e;
            live -= t.pipes[k].len;
        }
    }
    for (size_t k = 1; k < t.pipes.size(); k++) {
        if (t.pipes[k].rem == -1) { t.pipes[k].rem = t.m; t.aliveSum += t.pipes[k].len; }
    }
    return true;
}

static ll ceil20(ll x) { return (x + 19) / 20 * 20; }

// Независимая проверка раскладки: нет пересечений в пространстве и времени. ext = максимальный правый край.
static bool verify(const Test& t, const vector<int>& pos, ll& ext) {
    map<int, int> live;
    ext = 0;
    for (int e = 0; e < t.m; e++) {
        int c = t.cmds[e];
        if (c > 0) {
            int pid = t.evPid[e], p = pos[pid], L = c;
            if (p < 0) return false;
            auto it = live.lower_bound(p);
            if (it != live.end() && it->first < p + L) return false;
            if (it != live.begin()) { auto pr = prev(it); if (pr->first + pr->second > p) return false; }
            live[p] = L;
            ext = max(ext, (ll)p + L);
        } else {
            live.erase(pos[-c]);
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Эвристика: аллокатор по дыркам. oracle=false: обычный best-fit; oracle=true: знает сроки жизни.
// ---------------------------------------------------------------------------
struct Engine {
    const Test& t;
    bool oracle;
    double w2;
    int K;
    unordered_map<int, int> byStart, byEnd;   // дырки: start->end, end->start
    set<int> buckets[MAXL + 2];               // дырки по размеру, 41 = больше 40
    unordered_map<int, int> pAtStart, pAtEnd; // занятые участки: start->pid, end->pid
    vector<int> pos;
    ll trench = 0;
    bool overflow = false;

    Engine(const Test& tt, bool o, double w, int k) : t(tt), oracle(o), w2(w), K(k) {
        pos.assign(t.pipes.size(), -1);
    }
    static int bucketOf(int s) { return min(s, MAXL + 1); }
    void addFree(int s, int e) { byStart[s] = e; byEnd[e] = s; buckets[bucketOf(e - s)].insert(s); }
    int removeFree(int s) {
        auto it = byStart.find(s);
        if (it == byStart.end()) return s;
        int e = it->second;
        byStart.erase(it); byEnd.erase(e);
        buckets[bucketOf(e - s)].erase(s);
        return e;
    }

    void place(int pid) {
        int L = t.pipes[pid].len, remT = t.pipes[pid].rem;
        int bs = -1, bside = 0;
        double bscore = 1e300; int bD = 0;

        auto eval = [&](int s) {
            int e = byStart[s], h = e - s;
            for (int side = 0; side < 2; side++) {
                if (side == 1 && (h == L || e == (int)trench)) continue;   // у хвоста траншеи кладём только влево
                int D = side ? e - L : s;
                int nb = -1;
                if (side == 0) { if (s > 0) { auto it = pAtEnd.find(s); if (it != pAtEnd.end()) nb = it->second; } }
                else { auto it = pAtStart.find(e); if (it != pAtStart.end()) nb = it->second; }
                double mis = nb >= 0 ? fabs((double)t.pipes[nb].rem - remT) : 0.0;
                double score = (h - L) + w2 * mis / max(1, t.m);
                if (score < bscore - 1e-12 || (fabs(score - bscore) <= 1e-12 && D < bD)) {
                    bscore = score; bs = s; bside = side; bD = D;
                }
            }
        };

        if (!oracle) {
            for (int b = L; b <= MAXL + 1; b++)
                if (!buckets[b].empty()) { bs = *buckets[b].begin(); bside = 0; bD = bs; break; }
        } else {
            int seen = 0;
            for (int b = L; b <= MAXL + 1 && seen < K; b++)
                for (int s : buckets[b]) { eval(s); if (++seen >= K) break; }
        }

        int D;
        if (bs >= 0) {
            int e = removeFree(bs);
            D = bD;
            if (bside == 0) { if (bs + L < e) addFree(bs + L, e); }
            else { if (bs < D) addFree(bs, D); }
        } else {
            ll ts = trench;
            auto it = byEnd.find((int)trench);
            if (it != byEnd.end()) { ts = it->second; removeFree((int)ts); }
            ll need = L - (trench - ts);
            trench += 20 * ((need + 19) / 20);
            if (trench > MAX_TRENCH) { overflow = true; return; }
            D = (int)ts;
            if (ts + L < trench) addFree((int)(ts + L), (int)trench);
        }
        pos[pid] = D;
        pAtStart[D] = pid; pAtEnd[D + L] = pid;
    }

    void release(int pid) {
        int L = t.pipes[pid].len, s = pos[pid];
        pAtStart.erase(s); pAtEnd.erase(s + L);
        int ns = s, ne = s + L;
        auto li = byEnd.find(s);
        if (li != byEnd.end()) { int ls = li->second; removeFree(ls); ns = ls; }
        if (byStart.count(ne)) ne = removeFree(ne);
        addFree(ns, ne);
    }

    bool run() {
        for (int e = 0; e < t.m; e++) {
            int c = t.cmds[e];
            if (c > 0) { place(t.evPid[e]); if (overflow) return false; }
            else release(-c);
        }
        return true;
    }
};

// ---------------------------------------------------------------------------
// Точный поиск: можно ли уложить всё в траншею длины T (полный перебор позиций, события по порядку).
// ---------------------------------------------------------------------------
struct Exact {
    const Test& t;
    int T;
    ll nodes = 0, limit;
    bool aborted = false;
    vector<char> occ;
    vector<int> pos;

    Exact(const Test& tt, int T_, ll lim) : t(tt), T(T_), limit(lim) {
        occ.assign(T, 0);
        pos.assign(t.pipes.size(), -1);
    }

    bool dfs(int e) {
        if (e == t.m) return true;
        if (++nodes > limit) { aborted = true; return false; }
        int c = t.cmds[e];
        if (c < 0) {
            int k = -c, L = t.pipes[k].len, p = pos[k];
            for (int i = 0; i < L; i++) occ[p + i] = 0;
            if (dfs(e + 1)) return true;
            for (int i = 0; i < L; i++) occ[p + i] = 1;
            return false;
        }
        int L = c, pid = t.evPid[e];
        for (int p = 0; p + L <= T;) {
            int q = p + L - 1;
            while (q >= p && !occ[q]) q--;
            if (q >= p) { p = q + 1; continue; }     // окно занято: прыгаем за последнюю занятую клетку
            for (int i = 0; i < L; i++) occ[p + i] = 1;
            pos[pid] = p;
            if (dfs(e + 1)) return true;
            for (int i = 0; i < L; i++) occ[p + i] = 0;
            pos[pid] = -1;
            if (aborted) return false;
            p++;
        }
        return false;
    }
};

struct Result {
    ll T = 0;
    vector<int> pos;
    bool ok = false;
};

static Result heuristic(const Test& t, bool oracle, double w2, int K) {
    Engine en(t, oracle, w2, K);
    Result r;
    if (!en.run()) return r;
    ll ext;
    if (!verify(t, en.pos, ext)) return r;
    r.T = ceil20(ext);
    r.pos = en.pos;
    r.ok = true;
    return r;
}


// ---------------------------------------------------------------------------
// Эвристика зон: траншея фиксированной длины T. Долгоживущие трубы (срок жизни >= thr событий, либо
// "не удаляются никогда" при thr = LLONG_MAX) берут самую левую подходящую дырку и прижимаются влево;
// короткоживущие берут самую маленькую подходящую дырку (из равных самую правую) и прижимаются вправо.
// Так постоянные трубы копятся снизу, а временные занимают и освобождают верхнюю часть.
// ---------------------------------------------------------------------------
struct ZoneEngine {
    const Test& t;
    int T;
    ll thr;
    unordered_map<int, int> byStart, byEnd;
    set<int> buckets[MAXL + 2];
    vector<int> pos;

    ZoneEngine(const Test& tt, int T_, ll th) : t(tt), T(T_), thr(th) { pos.assign(t.pipes.size(), -1); }
    static int bucketOf(int s) { return min(s, MAXL + 1); }
    void addFree(int s, int e) { byStart[s] = e; byEnd[e] = s; buckets[bucketOf(e - s)].insert(s); }
    int removeFree(int s) {
        auto it = byStart.find(s);
        if (it == byStart.end()) return s;
        int e = it->second;
        byStart.erase(it); byEnd.erase(e);
        buckets[bucketOf(e - s)].erase(s);
        return e;
    }
    bool longLived(int pid) const {
        const Pipe& pp = t.pipes[pid];
        if (thr == LLONG_MAX) return pp.rem == t.m;
        return (ll)pp.rem - pp.laid >= thr;
    }
    bool place(int pid) {
        int L = t.pipes[pid].len;
        int s = -1, D = 0;
        if (longLived(pid)) {
            for (int b = L; b <= MAXL + 1; b++)
                if (!buckets[b].empty() && (s < 0 || *buckets[b].begin() < s)) s = *buckets[b].begin();
            if (s < 0) return false;
            int e = removeFree(s);
            D = s;
            if (s + L < e) addFree(s + L, e);
        } else {
            for (int b = L; b <= MAXL + 1; b++)
                if (!buckets[b].empty()) { s = *buckets[b].rbegin(); break; }
            if (s < 0) return false;
            int e = removeFree(s);
            D = e - L;
            if (s < D) addFree(s, D);
        }
        pos[pid] = D;
        return true;
    }
    void release(int pid) {
        int L = t.pipes[pid].len, s = pos[pid];
        int ns = s, ne = s + L;
        auto li = byEnd.find(s);
        if (li != byEnd.end()) { int ls = li->second; removeFree(ls); ns = ls; }
        if (byStart.count(ne)) ne = removeFree(ne);
        addFree(ns, ne);
    }
    bool run() {
        addFree(0, T);
        for (int e = 0; e < t.m; e++) {
            int c = t.cmds[e];
            if (c > 0) { if (!place(t.evPid[e])) return false; }
            else release(-c);
        }
        return true;
    }
};

// Ищет наименьшее T (кратное 20) из [LB, hiT - 20], при котором эвристика зон справляется. Монотонность по T предполагается.
static Result zoneSearch(const Test& t, ll LB, ll hiT, ll thr) {
    Result r;
    if (hiT - 20 < LB || LB <= 0) return r;
    auto feasible = [&](ll T, vector<int>* out) {
        ZoneEngine z(t, (int)T, thr);
        if (!z.run()) return false;
        if (out) *out = z.pos;
        return true;
    };
    ll lo = 0, hi = (hiT - 20 - LB) / 20;           // индексы k: T = LB + 20k
    vector<int> tmp;
    if (!feasible(LB + 20 * hi, &tmp)) return r;    // даже самая большая T из допустимых не помогает
    vector<int> bestPos = tmp;
    while (lo < hi) {
        ll mid = (lo + hi) / 2;
        if (feasible(LB + 20 * mid, &tmp)) { hi = mid; bestPos = tmp; }
        else lo = mid + 1;
    }
    ll ext;
    if (!verify(t, bestPos, ext)) return r;
    r.T = ceil20(ext); r.pos = bestPos; r.ok = true;
    return r;
}

int main(int argc, char** argv) {
    string dumpDir;
    double w2only = -1;
    int K = 48;
    int maxExactEvents = 400;
    ll nodeLimit = 20000000LL;
    vector<const char*> files;

    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "--dump" && i + 1 < argc) dumpDir = argv[++i];
        else if (a == "--w2" && i + 1 < argc) w2only = atof(argv[++i]);
        else if (a == "--k" && i + 1 < argc) K = atoi(argv[++i]);
        else if (a == "--max-exact-events" && i + 1 < argc) maxExactEvents = atoi(argv[++i]);
        else if (a == "--node-limit" && i + 1 < argc) nodeLimit = atoll(argv[++i]);
        else files.push_back(argv[i]);
    }
    if (files.empty()) {
        fprintf(stderr, "usage: %s [--dump dir] [--w2 X] [--k N] [--max-exact-events N] [--node-limit N] tests/*.in\n", argv[0]);
        return 1;
    }

    printf("%-8s %5s %7s %7s %8s %8s | %8s %7s | %8s %7s | %8s %7s %s\n",
           "test", "w", "pipes", "alive", "peak", "LB", "T_onl", "KPI", "T_orc", "KPI", "T_best", "KPI", "status");
    fflush(stdout);

    double wsum = 0, sUB = 0, sOn = 0, sBest = 0;

    for (const char* path : files) {
        Test t;
        if (!loadTest(path, t)) {
            fprintf(stderr, "[%s] %s\n", t.name.c_str(), t.err.c_str());
            continue;
        }
        if (t.hdrCount >= 0 && t.hdrCount != t.m)
            fprintf(stderr, "[%s] в заголовке %lld команд, прочитано %d\n", t.name.c_str(), t.hdrCount, t.m);

        int npipes = (int)t.pipes.size() - 1;
        ll LB = ceil20(t.peak);
        double kpiUB = LB > 0 ? (double)t.aliveSum / LB : 0;

        Result on = heuristic(t, false, 0, K);
        vector<double> w2s = w2only >= 0 ? vector<double>{w2only} : vector<double>{0, 10, 40, 160};
        Result orc;
        for (double w : w2s) {
            Result r = heuristic(t, true, w, K);
            if (r.ok && (!orc.ok || r.T < orc.T)) orc = r;
        }

        Result best = orc.ok ? orc : on;
        if (on.ok && (!best.ok || on.T < best.T)) best = on;

        // зоны: перебираем пороги долгоживучести, берём лучшее
        Result zon;
        string zonName;
        if (best.ok && best.T > LB) {
            vector<ll> thrs = {LLONG_MAX, 5000, 1000, 300, 100, 30};
            for (ll th : thrs) {
                Result r = zoneSearch(t, LB, best.T, th);
                if (r.ok && (!zon.ok || r.T < zon.T)) {
                    zon = r;
                    zonName = th == LLONG_MAX ? "навсегда" : ">=" + to_string(th);
                }
            }
            if (zon.ok && zon.T < best.T) best = zon; else zonName.clear();
        }

        string status, zpre;
        if (!zonName.empty() && zon.ok && best.T == zon.T) zpre = "[зоны " + zonName + "] ";
        status = zpre;
        if (!best.ok) status = "ОШИБКА эвристики";
        else if (best.T == LB) status += "OPT (T = LB)";
        else {
            status += "эвристика, зазор до LB " + to_string(best.T - LB) + " м";
            if (t.m <= maxExactEvents) {
                bool unknown = false, found = false;
                for (ll T = LB; T < best.T; T += 20) {
                    Exact ex(t, (int)T, nodeLimit);
                    bool ok = ex.dfs(0);
                    if (ok) {
                        ll ext;
                        if (verify(t, ex.pos, ext)) { best.T = ceil20(ext); best.pos = ex.pos; found = true; }
                        break;
                    }
                    if (ex.aborted) unknown = true;
                }
                if (found) status = zpre + (unknown ? "лучше эвристики, оптимум не доказан" : "OPT (перебор)");
                else status = zpre + (unknown ? "эвристика, перебор не завершён" : "OPT (перебор, эвристика оптимальна)");
            } else {
                status += "; перебор пропущен (много событий)";
            }
        }

        auto kpi = [&](const Result& r) { return r.ok && r.T > 0 ? (double)t.aliveSum / r.T : 0.0; };
        printf("%-8s %5g %7d %7lld %8lld %8lld | %8lld %7.4f | %8lld %7.4f | %8lld %7.4f %s\n",
               t.name.c_str(), t.weight, npipes, t.aliveSum, t.peak, LB,
               on.ok ? on.T : 0, kpi(on), orc.ok ? orc.T : 0, kpi(orc), best.ok ? best.T : 0, kpi(best), status.c_str());
        fflush(stdout);

        wsum += t.weight;
        sUB += t.weight * kpiUB;
        sOn += t.weight * kpi(on);
        sBest += t.weight * kpi(best);

        if (!dumpDir.empty() && best.ok) {
            string out = dumpDir + "/" + t.name + ".csv";
            FILE* f = fopen(out.c_str(), "w");
            if (!f) fprintf(stderr, "не удалось записать %s\n", out.c_str());
            else {
                fprintf(f, "pid,len,laid_event,removed_event,pos\n");
                for (int k = 1; k <= npipes; k++)
                    fprintf(f, "%d,%d,%d,%d,%d\n", k, t.pipes[k].len, t.pipes[k].laid,
                            t.pipes[k].rem == t.m ? -1 : t.pipes[k].rem, best.pos[k]);
                fclose(f);
            }
        }
    }

    printf("\nвзвешенная сумма KPI*w:  верхняя оценка (по LB) %.4f | online %.4f | лучшее найденное %.4f   (сумма весов %.2f)\n",
           sUB, sOn, sBest, wsum);
    return 0;
}