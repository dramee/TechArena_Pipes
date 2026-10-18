// Gas pipeline construction. Онлайн-решение.
//
// KPI = сумма длин оставшихся труб / длина траншеи. Числитель задан сценарием,
// значит задача = минимизировать выкопанное. Это динамическое распределение памяти.
//
// Режим 1 (обычный): best-fit с приоритетом точного совпадения по размеру.
//   Порог KMIN: предпочитаем дыру, оставляющую остаток >= KMIN, чтобы не плодить
//   огрызки по 1-3 м. На коротких тестах лучше KMIN=4, на длинных KMIN=16,
//   поэтому после 500-й трубы порог переключается (KMIN1=24: остаток >= 20 м ещё вмещает
//   типичную длинную трубу, остаток короче - мёртвая дыра).
//
// Режим 3 (предсказание периодичности): если длины входящих труб после последнего удаления образуют
//   периодическую последовательность (минимум PRIOR_REPS периодов подряд), считаем, что в таком потоке
//   удалят позиции с чётным номером внутри периода (так устроены открытые адверсарные тесты). Предполагаемые
//   к удалению трубы кладём вплотную к правому краю, остальные к левому, чтобы после удаления освободились
//   большие непрерывные блоки. Предположение проверяется по реальным удалениям и отключается навсегда,
//   если оно оказывается неверным чаще, чем верным.
//
// Сборка с -DSIM читает тест из stdin без интерактора и печатает итоговый KPI в stderr.
//
// Режим 2 (фазовый): включается, если обнаружена пачка из >= MIN_RUN труб одной
//   длины, в которой затем удалили часть труб по периодическому шаблону: каждый
//   остаток позиции по модулю p (p = 2..MAX_PERIOD) удалён целиком либо не тронут
//   (каждая вторая, каждая третья, две из трёх и т. д.). Тогда в следующих пачках:
//     - трубы, которые по выученному шаблону будут удалены, кладём к концу траншеи;
//     - остальные - к началу;
//     - копаем сразу на ожидаемый остаток пачки (среднее длин двух последних пачек), но
//       только если прогноз уже сбывался (в пределах 10%); иначе - не больше уложенного в пачке
//       (иначе после короткой пачки раскопанный про запас хвост остаётся пустым навсегда).
//   Так удаления освобождают один непрерывный блок вместо россыпи дыр размера k,
//   в которые трубы следующей пачки (размер k+1) физически не влезают.
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <utility>

#if defined(_WIN32) || defined(_MSC_VER)
  #define GETC getchar
#else
  #define GETC getchar_unlocked
#endif
using namespace std;
typedef long long ll;

static const int SM = 40;          // максимальная длина трубы
static const int MIN_RUN = 20;     // пачка короче не считается фазовой
static const int MAX_PERIOD = 8;
#ifndef PRIOR_PERIODIC
#define PRIOR_PERIODIC 1   // предсказание: в периодической последовательности длин удаляют позиции с чётным номером внутри периода
#endif
static const int PRIOR_MAX_PERIOD = 12;
static const int PRIOR_MIN_PERIOD = 3;   // меньшие периоды слишком часто возникают случайно
static const int PRIOR_REPS = 6;         // сколько периодов подряд должно совпасть, чтобы поверить в периодичность   // максимальный период шаблона удалений внутри серии
#ifndef SWITCH_AT_V
#define SWITCH_AT_V 500
#endif
#ifndef KMIN0
#define KMIN0 4
#endif
#ifndef KMIN1
#define KMIN1 24
#endif
#ifndef CAP_MULT
#define CAP_MULT 1
#endif
#ifndef GUARD_V
#define GUARD_V 100
#endif
static const ll SWITCH_AT = SWITCH_AT_V;   // после скольких труб поднять порог остатка
static ll KMIN = KMIN0;
static const ll GUARD = GUARD_V;   // правило правой дыры работает только на коротких траншеях
static ll trenchLen = 0;       // текущая длина траншеи (нужна внутри goodFit)

static map<ll,ll> fa, fe;          // дыры: start -> len, end -> start
static set<pair<ll,ll>> bySize;    // (len, start)
static set<ll> ba[SM+2];           // дыры по размеру: 1..40, корзина 41 = больше 40

static inline int bk(ll l) { return (int)(l < SM + 1 ? l : SM + 1); }

static inline void addRaw(ll s, ll len) {
    if (len <= 0) return;
    fa[s] = len; fe[s+len] = s;
    bySize.insert({len, s}); ba[bk(len)].insert(s);
}
static inline void delRaw(ll s, ll len) {
    fa.erase(s); fe.erase(s+len);
    bySize.erase({len, s}); ba[bk(len)].erase(s);
}
static void addFree(ll s, ll len) {           // освободить участок со склейкой соседей
    auto a = fe.find(s);
    if (a != fe.end()) { ll ps = a->second, pl = s - ps; delRaw(ps, pl); s = ps; len += pl; }
    auto b = fa.find(s + len);
    if (b != fa.end()) { ll nl = b->second; delRaw(s + len, nl); len += nl; }
    addRaw(s, len);
}

static bool rightFit(ll L, ll& hs, ll& hl);

// best-fit с порогом остатка, минуя дыру с началом ban
static bool pickExcept(ll L, ll ban, ll& hs, ll& hl) {   // best-fit с порогом, минуя дыру ban
    for (int pass = 0; pass < 2; pass++) {
        ll from = pass == 0 ? L + KMIN : L;
        for (auto it = bySize.lower_bound({from, LLONG_MIN}); it != bySize.end(); ++it) {
            if (it->second == ban) continue;
            hl = it->first; hs = it->second; return true;
        }
    }
    return false;
}
// Выбор дыры в обычном режиме.
// Точное совпадение по размеру — лучший вариант, но есть исключение: если точная дыра
// оказалась самой правой, её лучше поберечь. Дыры у конца траншеи чаще сливаются с
// соседними при удалениях, и занимать их, когда есть альтернатива, — значит терять
// будущий большой свободный участок.
// первая непустая корзина размеров [lo, hi] (по возрастанию); корзина 41 - любые дыры больше 40
static bool firstBucket(ll lo, ll hi, ll ban, ll& hs, ll& hl) {
    for (ll k = lo; k <= hi && k <= SM; k++) {
        if (ba[k].empty()) continue;
        if (ba[k].size() == 1 && *ba[k].begin() == ban) continue;
        ll st = *ba[k].begin();                          // из равных по размеру берём левую
        if (st == ban) st = *next(ba[k].begin());
        hs = st; hl = k; return true;
    }
    if (hi > SM) {
        for (auto it = bySize.lower_bound({SM + 1, LLONG_MIN}); it != bySize.end(); ++it) {
            if (it->second == ban) continue;
            hl = it->first; hs = it->second; return true;
        }
    }
    return false;
}
static bool goodFit(ll L, ll& hs, ll& hl) {
    ll rs, rl;
    bool anyRight = rightFit(L, rs, rl);                  // самая правая подходящая дыра
    ll ban = -1;
    if (!ba[(int)L].empty()) {                            // есть точное совпадение
        ll es = *ba[(int)L].begin();
        if (!anyRight || es != rs || trenchLen > GUARD) { hs = es; hl = L; return true; }
        ban = es;                                         // точная дыра самая правая - поищем другую
    }
    if (firstBucket(L + KMIN, LLONG_MAX, ban, hs, hl)) return true;
    if (firstBucket(L, L + KMIN - 1, ban, hs, hl)) return true;
    if (ban >= 0) { hs = ban; hl = L; return true; }
    return false;
}
static bool leftFit(ll L, ll& hs, ll& hl) {   // самая левая подходящая дыра
    ll b = LLONG_MAX; int bl = -1;
    for (int k = (int)L; k <= SM + 1; k++)
        if (!ba[k].empty() && *ba[k].begin() < b) { b = *ba[k].begin(); bl = k; }
    if (bl < 0) return false;
    hs = b; hl = fa[b]; return true;
}
static bool rightFit(ll L, ll& hs, ll& hl) {  // самая правая подходящая дыра
    ll b = -1; int bl = -1;
    for (int k = (int)L; k <= SM + 1; k++)
        if (!ba[k].empty() && *ba[k].rbegin() > b) { b = *ba[k].rbegin(); bl = k; }
    if (bl < 0) return false;
    hs = b; hl = fa[b]; return true;
}

// Наименьший период p (2..PRIOR_MAX_PERIOD), которому подчиняются последние 3p длин, при условии что в периоде есть разные длины.
static int detectPeriod(const vector<int>& v) {
    int n = (int)v.size();
    for (int p = PRIOR_MIN_PERIOD; p <= PRIOR_MAX_PERIOD; p++) {
        if (n < PRIOR_REPS * p) break;
        bool ok = true;
        for (int i = n - PRIOR_REPS * p; i + p < n && ok; i++) if (v[i] != v[i + p]) ok = false;
        if (!ok) continue;
        bool allSame = true;
        for (int i = n - p + 1; i < n; i++) if (v[i] != v[n - p]) allSame = false;
        if (!allSame) return p;
    }
    return 0;
}

static inline int readInt() {
    int c = GETC();
    while (c != '-' && (c < '0' || c > '9')) { if (c == EOF) return INT_MIN; c = GETC(); }
    int sg = 1; if (c == '-') { sg = -1; c = GETC(); }
    int v = 0; while (c >= '0' && c <= '9') { v = v * 10 + (c - '0'); c = GETC(); }
    return v * sg;
}

int main() {
    ll& trench = trenchLen;
    ll laid = 0;
#ifdef SIM   // локальный режим: тест читается из stdin без интерактора (заголовок "вес n" пропускаем)
    { double w; int n; if (scanf("%lf %d", &w, &n) != 2) return 1; }
#endif
    vector<pair<ll,ll>> pipes;  pipes.reserve(200005);  pipes.push_back({-1,-1});
    vector<int> runIdx;         runIdx.reserve(200005); runIdx.push_back(-1);
    vector<int> runId;          runId.reserve(200005);  runId.push_back(-1);

    int curSize = -1, curRun = 0, curCount = 0;
    vector<int> runLens;                       // длины укладок после последнего удаления (для поиска периода)
    vector<signed char> predDoom;              // предсказание судьбы трубы: -1 нет, 0 останется, 1 удалят
    predDoom.reserve(200005); predDoom.push_back(-1);
    ll priorHit = 0, priorMiss = 0; bool priorBad = false;
    vector<int> remIdx;                        // позиции удалённых труб текущей пачки
    bool havePattern = false, expTrusted = false;
    int patPeriod = 2, patMask = 0, patRunLen = 0, patRunLenPrev = 0, expLen = 0;   // шаблон: удаляются позиции i, у которых (patMask >> (i % patPeriod)) & 1

    for (;;) {
        int cmd = readInt();
        if (cmd == 0 || cmd == INT_MIN) break;

        if (cmd > 0) {
            ll L = cmd;

            if ((int)L != curSize) {           // пачка сменилась: анализируем предыдущую
                if (curSize >= 0 && curCount >= MIN_RUN && !remIdx.empty()) {
                    // прогноз длины пачки заслуживает доверия, если прошлый прогноз сбылся с точностью 10%
                    expTrusted = havePattern && expLen > 0 && llabs((ll)curCount - expLen) * 10 <= expLen;
                    ll cnt = (ll)remIdx.size();
                    bool found = false;
                    if (cnt * 4 >= curCount && cnt * 4 <= (ll)curCount * 3) {
                        // ищем наименьший период p, при котором каждый остаток по модулю p либо удалён целиком, либо не тронут
                        for (int p = 2; p <= MAX_PERIOD && !found; p++) {
                            int remC[MAX_PERIOD] = {0}, totC[MAX_PERIOD] = {0};
                            for (int i = 0; i < curCount; i++) totC[i % p]++;
                            for (size_t i = 0; i < remIdx.size(); i++) remC[remIdx[i] % p]++;
                            bool ok = true; int mask = 0;
                            for (int r = 0; r < p && ok; r++) {
                                if (remC[r] == totC[r] && totC[r] > 0) mask |= 1 << r;
                                else if (remC[r] != 0) ok = false;
                            }
                            if (ok && mask != 0 && mask != (1 << p) - 1) {
                                found = true; patPeriod = p; patMask = mask;
                                patRunLenPrev = havePattern ? patRunLen : curCount; patRunLen = curCount;
                                expLen = (patRunLen + patRunLenPrev) / 2;   // ожидаемая длина следующей серии (среднее двух последних)
                            }
                        }
                    }
                    havePattern = found;
                }
                curSize = (int)L; curRun++; curCount = 0; remIdx.clear();
            }

            if (++laid == SWITCH_AT) KMIN = KMIN1;
            int idx = curCount++;
#ifdef NOPHASE
            havePattern = false;
#endif
            bool phase  = havePattern && patRunLen >= MIN_RUN;
            bool doomed = phase && (((patMask >> (idx % patPeriod)) & 1) != 0);
            bool prior = false;                // решение принято по предположению, а не по выученному шаблону
            runLens.push_back((int)L);
#ifdef NOPHASE
            priorBad = true;
#endif
            if (!phase && !priorBad) {
#if PRIOR_PERIODIC
                int per = detectPeriod(runLens);
                if (per > 0) {
                    int res = (int)(runLens.size() - 1) % per;
                    phase = true; prior = true;
                    doomed = (res % 2 == 0) && !(per % 2 == 1 && res == per - 1);
                }
#endif
            }

            ll hs = -1, hl = -1, D, dig = 0;
            bool got = !phase ? goodFit(L, hs, hl)
                              : (doomed ? rightFit(L, hs, hl) : leftFit(L, hs, hl));
            if (got) {
                delRaw(hs, hl);
                if (!doomed) { D = hs;           if (hl > L) addRaw(hs + L, hl - L); }
                else         { D = hs + hl - L;   if (hl > L) addRaw(hs, hl - L); }
            } else {
                ll ts = trench, tl = 0;
                auto e = fe.find(trench);
                if (e != fe.end()) { ts = e->second; tl = trench - ts; delRaw(ts, tl); }
                ll need = L - tl;
                if (phase && !prior) {            // копаем сразу на остаток пачки
                    ll rest = (ll)(expLen - idx) * L - tl;
#if CAP_MULT > 0   // пока прогноз длины пачки не подтвердился, не копаем вперёд больше CAP_MULT x уже уложенного в пачке
                    ll cap = (ll)CAP_MULT * (ll)idx * L + 40;
                    if (!expTrusted && rest > cap) rest = cap;
#endif
                    if (rest > need) need = rest;
                }
                dig = (need + 19) / 20;
                trench += 20 * dig;
                ll nl = tl + 20 * dig;
                if (!doomed) { D = ts;             if (nl > L) addRaw(ts + L, nl - L); }
                else         { D = ts + nl - L;    if (nl > L) addRaw(ts, nl - L); }
            }

            pipes.push_back({D, L}); runIdx.push_back(idx); runId.push_back(curRun);
            predDoom.push_back(prior ? (doomed ? 1 : 0) : -1);
            printf("%lld %lld\n", dig, D);
            fflush(stdout);
        } else {
            ll k = -cmd;
            printf("%lld\n", pipes[k].first);
            fflush(stdout);
            addFree(pipes[k].first, pipes[k].second);
            if (runId[k] == curRun) remIdx.push_back(runIdx[k]);
            if (predDoom[k] >= 0) {                   // сверяем предположение с реальностью
                if (predDoom[k] == 1) priorHit++; else priorMiss++;
                if (priorHit + priorMiss >= 50 && priorMiss > priorHit) priorBad = true;
            }
            runLens.clear();
            pipes[k] = {-1, -1};
        }
    }
#ifdef SIM
    { ll alive = 0; for (size_t i = 1; i < pipes.size(); i++) if (pipes[i].first >= 0) alive += pipes[i].second;
      fprintf(stderr, "RATIO %.6f alive %lld trench %lld\n", (double)alive / (double)trench, alive, trench); }
#endif
    return 0;
}