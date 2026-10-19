// sim.cpp: локальный симулятор политик укладки труб (без интерактора).
// Сборка: clang++ -O2 -std=c++17 sim.cpp -o sim     (или g++)
// Запуск: ./sim tests/1.in tests/2.in ...           (или ./sim tests/*.in)
//         ./sim --only best-left,s3:BB tests/*.in   (только выбранные политики)
// Отладка: clang++ -g -O1 -std=c++17 -fsanitize=address,undefined sim.cpp -o sim_dbg
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

enum Mode { BEST, FIRST, WORST, SMALL_BIG, PCT_BIG, LAST, SIGMA3, BESTK, MT };

struct Cfg {
    string name;
    Mode mode;
    bool right;          // прижимать трубу к правому краю дырки
    double eps;          // SMALL_BIG: "маленькая", если L < eps * stat
    long long switchAt;  // SMALL_BIG: до этого числа труб среднее, дальше медиана
    double pct;          // PCT_BIG: доля прежних труб короче L, с которой берём самую большую дырку
    // SIGMA3: трубы делятся на классы по правилу 3 сигм. Интервал [mu-3s, mu+3s] режется на три
    // равные части шириной 2s, т.е. границы классов mu - kSigma*s и mu + kSigma*s (kSigma = 1).
    // Маленькие -> worst-fit, средние -> medMode, большие -> bigMode.
    Mode medMode = BEST;
    Mode bigMode = BEST;
    double kSigma = 1.0;
    double alpha = 0.0;  // 0 - статистика по всем трубам; иначе экспоненциальное сглаживание с весом alpha
    bool smallNoTail = false;  // worst-fit для маленьких без хвостовой дырки у конца траншеи
    // Условие включения worst-fit для маленьких (иначе маленькая труба идёт по политике средних):
    //   GATE_NONE  - всегда;
    //   GATE_GAP   - между самой длинной маленькой и самой короткой не маленькой длиной есть
    //                хотя бы gapMin ни разу не встреченных длин (маленькие - отдельный кластер);
    //   GATE_RARE  - маленькие встречаются не чаще rareShare и состоят не более чем из rareVals длин;
    //   GATE_FEW   - во всём потоке не больше fewVals разных длин (дискретный, «адверсарный» поток);
    //   GATE_SHORT - пока уложено меньше shortN труб.
    int gate = 0;
    int gapMin = 1;
    double rareShare = 0.35;
    int rareVals = 1;
    int fewVals = 6;
    long long shortN = 100;
    // MT: классы Мартелло-Тота относительно дырки h (вместимость C = h) с адаптивным alpha:
    //   alpha = квантиль mtQ (в %) длин последних mtWin труб - короче alpha трубы почти не приходят;
    //   остаток r = h - L: r = 0 - точное заполнение; r < alpha - труба «большая» для дырки, остаток мёртвый;
    //   r >= alpha - остаток ещё пригоден. Пары (по процедуре сведения MTRP): остаток r сам является
    //   частой длиной (доля r среди недавних >= mtTau), т.е. дырку скорее всего закроют две трубы без хвоста.
    //   Приоритет: точное > [пара, если mtPair] > пригодный остаток (наименьшая дырка) > мёртвый (наименьший r).
    //   mtPair: 0 - без пар, 1 - пара выше пригодного остатка, 2 - пара только вместо мёртвого остатка.
    int mtQ = 10;
    double mtTau = 0.05;
    int mtPair = 1;
    int mtWin = 2000;
};
enum Gate { GATE_NONE, GATE_GAP, GATE_RARE, GATE_FEW, GATE_SHORT };

static const long long MAX_TRENCH = 100000000LL;
static const int MAXL = 40;

struct Allocator {
    Cfg cfg;
    unordered_map<int, int> byStart, byEnd;  // start -> end, end -> start
    set<pair<int, int>> bySize;              // (размер, start)
    set<int> buckets[MAXL + 2];              // размеры 1..40, корзина 41 = больше 40
    long long lenCount[MAXL + 1];
    long long cnt = 0, sumLen = 0, trench = 0;
    double mu = 0, m2 = 0;                   // SIGMA3: среднее и сумма квадратов отклонений (Уэлфорд)
    long long smallCnt = 0, medCnt = 0, bigCnt = 0;
    vector<int> win;                          // MT: последние mtWin длин
    long long winPos = 0, wCnt[MAXL + 2] = {0}, wTot = 0;

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
    int lastFit(int L) const {   // самая правая подходящая дырка
        int s = -1;
        for (int b = L; b <= MAXL + 1; b++)
            if (!buckets[b].empty() && *buckets[b].rbegin() > s) s = *buckets[b].rbegin();
        return s;
    }
    // best-fit с порогом остатка: точное совпадение, иначе дырка с остатком >= kmin, иначе любая
    // (kmin = 4 первые 500 труб, дальше 24 - как в solution.cpp)
    int bestFitK(int L) const {
        if (!buckets[L].empty()) return *buckets[L].begin();
        int kmin = cnt < 500 ? 4 : 24;
        auto it = bySize.lower_bound(make_pair(L + kmin, -1));
        if (it != bySize.end()) return it->second;
        return bestFit(L);
    }
    // самая большая дырка, кроме хвостовой (той, что упирается в конец траншеи)
    int worstFitNoTail(int L) const {
        for (auto it = bySize.rbegin(); it != bySize.rend() && it->first >= L; ++it) {
            if (it->second + it->first == trench) continue;
            return it->second;
        }
        return worstFit(L);
    }
    int mtAlpha() const {
        if (wTot == 0) return 1;
        long long need = (wTot * cfg.mtQ + 99) / 100, acc = 0;
        for (int v = 1; v <= MAXL; v++) { acc += wCnt[v]; if (acc >= need) return v; }
        return MAXL;
    }
    int mtFit(int L) const {
        if (!buckets[L].empty()) return *buckets[L].begin();          // точное заполнение
        int alpha = mtAlpha();
        int pair = -1, reuse = -1, dead = -1;
        double pairFreq = -1;
        for (int k = L + 1; k <= MAXL; k++) {
            if (buckets[k].empty()) continue;
            int r = k - L, st = *buckets[k].begin();
            double fr = wTot ? (double)wCnt[r] / (double)wTot : 0;
            if (cfg.mtPair && fr >= cfg.mtTau && (cfg.mtPair == 1 || r < alpha) && fr > pairFreq) { pairFreq = fr; pair = st; }
            if (r >= alpha) { if (reuse < 0) reuse = st; }
            else if (dead < 0) dead = st;
        }
        if (reuse < 0) {                                              // дырки длиннее 40 - остаток всегда пригоден
            auto it = bySize.lower_bound(make_pair(MAXL + 1, -1));
            if (it != bySize.end()) reuse = it->second;
        }
        if (cfg.mtPair == 1 && pair >= 0) return pair;
        if (reuse >= 0) return reuse;
        if (pair >= 0) return pair;
        return dead;
    }
    int byMode(Mode m, int L) const {
        switch (m) {
            case BESTK: return bestFitK(L);
            case FIRST: return firstFit(L);
            case WORST: return worstFit(L);
            case LAST: return lastFit(L);
            default: return bestFit(L);
        }
    }
    // приближённое СКО по уже уложенным трубам
    double sigma() const {
        if (cfg.alpha > 0) return sqrt(max(m2, 0.0));
        return cnt > 1 ? sqrt(m2 / (double)(cnt - 1)) : 0.0;
    }
    // 0 - маленькая, 1 - средняя, 2 - большая (пока статистики мало - средняя)
    int sizeClass(int L) const {
        if (cnt < 2) return 1;
        double s = sigma(), lo = mu - cfg.kSigma * s, hi = mu + cfg.kSigma * s;
        if (s <= 0) return 1;
        return L < lo ? 0 : (L > hi ? 2 : 1);
    }
    bool gateOpen() const {
        if (cfg.gate == GATE_NONE) return true;
        if (cfg.gate == GATE_SHORT) return cnt < cfg.shortN;
        double lo = mu - cfg.kSigma * sigma();
        int distinct = 0, smallVals = 0, maxSmall = 0, minRest = MAXL + 1;
        long long smallCnt = 0;
        for (int v = 1; v <= MAXL; v++) {
            if (!lenCount[v]) continue;
            distinct++;
            if (v < lo) { smallVals++; smallCnt += lenCount[v]; maxSmall = v; }
            else if (minRest > MAXL) minRest = v;
        }
        switch (cfg.gate) {
            case GATE_GAP: return maxSmall > 0 && minRest <= MAXL && minRest - maxSmall - 1 >= cfg.gapMin;
            case GATE_RARE: return smallVals > 0 && smallVals <= cfg.rareVals && (double)smallCnt <= cfg.rareShare * (double)cnt;
            case GATE_FEW: return distinct <= cfg.fewVals;
        }
        return true;
    }
    void updateStats(int L) {
        if (cfg.alpha > 0) {      // экспоненциально сглаженные среднее и дисперсия
            if (cnt == 0) { mu = L; m2 = 0; return; }
            double d = L - mu;
            mu += cfg.alpha * d;
            m2 = (1 - cfg.alpha) * (m2 + cfg.alpha * d * d);
        } else {                  // алгоритм Уэлфорда (cnt уже учитывает новую трубу)
            double d = L - mu;
            mu += d / (double)cnt;
            m2 += d * (L - mu);
        }
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
            case LAST: return lastFit(L);
            case BESTK: return bestFitK(L);
            case MT: return mtFit(L);
            case SIGMA3: {
                int c = sizeClass(L);
                if (c == 0 && gateOpen()) return cfg.smallNoTail ? worstFitNoTail(L) : worstFit(L);
                if (c == 0) return byMode(cfg.medMode, L);
                return byMode(c == 1 ? cfg.medMode : cfg.bigMode, L);
            }
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
        if (cfg.mode == SIGMA3) { int c = sizeClass(L); (c == 0 ? smallCnt : c == 1 ? medCnt : bigCnt)++; }
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
        updateStats(L);
        if (cfg.mode == MT) {
            if ((int)win.size() < cfg.mtWin) win.push_back(L);
            else { wCnt[win[winPos % cfg.mtWin]]--; wTot--; win[winPos % cfg.mtWin] = L; }
            winPos++; wCnt[L]++; wTot++;
        }
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
    // SIGMA3: маленькие всегда worst-fit, для средних и больших перебираем политики
    const pair<Mode, const char*> pol[] = {{BEST, "B"}, {BESTK, "K"}, {FIRST, "F"}, {WORST, "W"}, {LAST, "L"}};
    for (auto& m : pol)
        for (auto& b : pol) {
            Cfg c{string("s3:") + m.second + b.second, SIGMA3, false, 0, 0, 0};
            c.medMode = m.first; c.bigMode = b.first;
            cfgs.push_back(c);
        }
    for (auto mb : {make_pair(BEST, BEST), make_pair(BESTK, BESTK), make_pair(BEST, BESTK), make_pair(BESTK, BEST)})
        for (double k : {0.5, 1.25, 1.5, 1.75, 2.0}) {
            char nm[32]; snprintf(nm, sizeof nm, "s3:%c%c k%.2f", mb.first == BEST ? 'B' : 'K', mb.second == BEST ? 'B' : 'K', k);
            Cfg c{nm, SIGMA3, false, 0, 0, 0}; c.medMode = mb.first; c.bigMode = mb.second; c.kSigma = k; cfgs.push_back(c);
        }
    for (double a : {0.01, 0.001}) {
        char nm[32]; snprintf(nm, sizeof nm, "s3:KK a%g", a);
        Cfg c{nm, SIGMA3, false, 0, 0, 0}; c.medMode = c.bigMode = BESTK; c.alpha = a; cfgs.push_back(c);
    }
    {   // вариант: worst-fit для маленьких берёт самую большую дырку, но не хвост траншеи
        Cfg c{"s3:KK notail", SIGMA3, false, 0, 0, 0}; c.medMode = c.bigMode = BESTK; c.smallNoTail = true; cfgs.push_back(c);
        Cfg d{"s3:BB notail", SIGMA3, false, 0, 0, 0}; d.smallNoTail = true; cfgs.push_back(d);
    }
    cfgs.push_back({"bestK-left", BESTK, false, 0, 0, 0});
    // SIGMA3 KK, k=1 с разными условиями включения worst-fit для маленьких
    auto gated = [&](const char* nm, int gate, auto tweak) {
        Cfg c{nm, SIGMA3, false, 0, 0, 0}; c.medMode = c.bigMode = BESTK; c.gate = gate; tweak(c); cfgs.push_back(c);
    };
    gated("g:gap1", GATE_GAP, [](Cfg& c) { c.gapMin = 1; });
    gated("g:gap2", GATE_GAP, [](Cfg& c) { c.gapMin = 2; });
    gated("g:gap3", GATE_GAP, [](Cfg& c) { c.gapMin = 3; });
    gated("g:rare1/35", GATE_RARE, [](Cfg& c) { c.rareVals = 1; c.rareShare = 0.35; });
    gated("g:rare1/20", GATE_RARE, [](Cfg& c) { c.rareVals = 1; c.rareShare = 0.20; });
    gated("g:rare2/35", GATE_RARE, [](Cfg& c) { c.rareVals = 2; c.rareShare = 0.35; });
    gated("g:few4", GATE_FEW, [](Cfg& c) { c.fewVals = 4; });
    gated("g:few6", GATE_FEW, [](Cfg& c) { c.fewVals = 6; });
    gated("g:few8", GATE_FEW, [](Cfg& c) { c.fewVals = 8; });
    gated("g:short50", GATE_SHORT, [](Cfg& c) { c.shortN = 50; });
    gated("g:short200", GATE_SHORT, [](Cfg& c) { c.shortN = 200; });
    // Мартелло-Тот: классы относительно дырки + пары MTRP
    auto mt = [&](int q, int pairMode, double tau, int w) {
        char nm[40];
        snprintf(nm, sizeof nm, "mt:q%d%s%s", q, pairMode == 0 ? "" : pairMode == 1 ? "p" : "d",
                 pairMode ? (tau < 0.04 ? "2" : tau < 0.08 ? "5" : "10") : "");
        string name = nm; if (w != 2000) name += "w" + to_string(w);
        Cfg c{name, MT, false, 0, 0, 0}; c.mtQ = q; c.mtPair = pairMode; c.mtTau = tau; c.mtWin = w; cfgs.push_back(c);
    };
    for (int q : {1, 5, 10, 20, 35}) mt(q, 0, 0, 2000);
    for (int q : {5, 10, 20})
        for (int pm : {1, 2})
            for (double tau : {0.02, 0.05, 0.10}) mt(q, pm, tau, 2000);
    mt(10, 0, 0, 300); mt(10, 0, 0, 20000); mt(10, 2, 0.05, 300); mt(10, 2, 0.05, 20000);

    // --only name1,name2,... - запустить только перечисленные политики
    int argBase = 1;
    if (argc > 2 && strcmp(argv[1], "--only") == 0) {
        string list = string(",") + argv[2] + ",";
        vector<Cfg> keep;
        for (auto& c : cfgs) if (list.find("," + c.name + ",") != string::npos) keep.push_back(c);
        cfgs = keep;
        argBase = 3;
    }

    if (argc <= argBase) {
        fprintf(stderr, "usage: %s [--only name,name] tests/*.in\n", argv[0]);
        return 1;
    }

    printf("%-22s", "test");
    for (auto& c : cfgs) printf("%-14s  ", c.name.c_str());
    printf("\n");
    fflush(stdout);

    vector<double> total(cfgs.size(), 0.0);   // сумма KPI
    vector<double> wtotal(cfgs.size(), 0.0);  // сумма KPI * вес
    double wsum = 0.0;                        // сумма весов
    int files = 0;

    for (int i = argBase; i < argc; i++) {
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
                printf("%-14s  ", "INVALID");
                fprintf(stderr, "\n[%s / %s] %s\n", name.c_str(), cfgs[j].name.c_str(), why.c_str());
            } else {
                printf("%-14.4f  ", kpi);
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
    for (size_t j = 0; j < cfgs.size(); j++) printf("%-14.4f  ", files ? total[j] / files : 0.0);
    printf("\n");
    printf("%-22s", "weighted mean");
    for (size_t j = 0; j < cfgs.size(); j++) printf("%-14.4f  ", wsum > 0 ? wtotal[j] / wsum : 0.0);
    printf("\n");
    printf("%-22s", "score (sum KPI*w)");
    for (size_t j = 0; j < cfgs.size(); j++) printf("%-14.4f  ", wtotal[j]);
    printf("\n");
    return 0;
}