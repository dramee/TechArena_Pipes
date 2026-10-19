// Gas pipeline construction. Онлайн-решение.
//
// KPI = сумма длин оставшихся труб / длина траншеи. Числитель задан сценарием,
// значит задача = минимизировать выкопанное. Это динамическое распределение памяти.
//
// Режим 1 (обычный): best-fit с приоритетом точного совпадения по размеру.
//   Порог KMIN: предпочитаем дыру, оставляющую остаток >= KMIN, чтобы не плодить
//   огрызки по 1-3 м. На коротких тестах лучше KMIN=4, на длинных KMIN=24 (остаток >= 20 м
//   ещё вмещает типичную длинную трубу, остаток короче - мёртвая дыра), поэтому после
//   500-й трубы порог переключается.
//
// Режим 2 (фазовый): включается, если обнаружена пачка из >= MIN_RUN труб одной
//   длины (пачка кончается при смене длины или после всплеска удалений), в которой затем
//   удалили часть труб по периодическому шаблону: каждый остаток позиции по модулю p
//   (p = 2..MAX_PERIOD) удалён целиком либо не тронут (каждая вторая, каждая третья,
//   две из трёх и т. д.). Тогда в следующих пачках:
//     - трубы, которые по выученному шаблону будут удалены, кладём к концу траншеи;
//     - остальные - к началу;
//     - копаем сразу на ожидаемый остаток пачки (среднее длин двух последних пачек), но
//       только если прогноз уже сбывался (в пределах 10%); иначе - не больше уложенного в пачке
//       (иначе после короткой пачки раскопанный про запас хвост остаётся пустым навсегда).
//   Так удаления освобождают один непрерывный блок вместо россыпи дыр размера k,
//   в которые трубы следующей пачки (размер k+1) физически не влезают.
//
// Режим 3 (предсказание периодичности): если длины входящих труб после последнего удаления образуют
//   периодическую последовательность (минимум PRIOR_REPS периодов подряд), считаем, что в таком потоке
//   удалят позиции с чётным номером внутри периода (так устроены открытые адверсарные тесты). Предполагаемые
//   к удалению трубы кладём вплотную к правому краю, остальные к левому, чтобы после удаления освободились
//   большие непрерывные блоки. Предположение проверяется по реальным удалениям и отключается навсегда,
//   если оно оказывается неверным чаще, чем верным.
//
// Режим 4 (выученная судьба позиций периода): период ищется по всем укладкам подряд
//   (удаления его не сбрасывают). Каждую трубу «судим» через L4_H укладок: удалили её за это время
//   или нет. Если каждая позиция периода почти всегда удаляется либо почти никогда (>= 80% / <= 20%),
//   кладём по выученному шаблону:
//     - остающиеся - в дыру точно по размеру, иначе в хвост траншеи (плотный рост без дыр);
//     - обречённые - best-fit по дырам: они переиспользуют места прошлых обречённых;
//       хвост копаем сразу под все обречённые трубы текущего периода.
//   Главное - остающиеся не занимают и не дробят дыры, которые освобождают обречённые.
//   Если прогнозы ошибаются чаще чем в 20% случаев, режим выключается навсегда.
//
// Сборка с -DSIM читает тест из stdin без интерактора и печатает итоговый KPI в stderr.
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <vector>
#include <set>
#include <utility>

#if defined(_WIN32)
  #define GETC _getchar_nolock
#else
  #define GETC getchar_unlocked
#endif
using namespace std;
typedef long long ll;

static const int SM = 40;          // максимальная длина трубы
static const int MIN_RUN = 20;     // пачка короче не считается фазовой
static const int MAX_PERIOD = 8;   // максимальный период шаблона удалений внутри пачки
static const int PRIOR_MAX_PERIOD = 12;
static const int PRIOR_MIN_PERIOD = 3;   // меньшие периоды слишком часто возникают случайно
static const int PRIOR_REPS = 6;         // сколько периодов подряд должно совпасть, чтобы поверить в периодичность
static const ll SWITCH_AT = 500;         // после скольких труб поднять порог остатка
static const ll KMIN1 = 24;
static const ll GUARD = 100;             // правило правой дыры работает только на коротких траншеях
static const ll CAP_MULT = 1;            // без подтверждённого прогноза копаем вперёд не больше CAP_MULT x уложенного в пачке
static const int L4_H = 64;              // режим 4: через сколько укладок трубу считаем оставшейся навсегда
static const int L4_PURE = 80;           // режим 4: позиция чистая, если её удаляют в >= 80% случаев или в <= 20%
static ll KMIN = 4;
static ll trenchLen = 0;                 // текущая длина траншеи (нужна внутри goodFit)

// Дыры. Позиции и длины укладываются в int: интерактор ограничивает траншею 1e8 м.
//   hLen(s) - длина дыры, начинающейся в s (0 - дыры нет);
//   hBeg(e) - начало+1 дыры, кончающейся в e (0 - дыры нет);
//   ba[k]   - начала дыр длины k (1..40), ba[41] - начала дыр длиннее 40;
//   big     - дыры длиннее 40, упорядоченные по (длина, начало).
// hLen/hBeg - массивы по позиции, разбитые на страницы по 64K: память растёт вместе с траншеей
// (8 байт на метр), без перевыделений и копирования при росте.
struct PagedArray {
    static const int SH = 16, PS = 1 << SH;
    vector<int*> pages;
    ~PagedArray() { for (int* p : pages) free(p); }
    int get(ll i) const { size_t k = (size_t)(i >> SH); return k < pages.size() && pages[k] ? pages[k][i & (PS - 1)] : 0; }
    void set(ll i, int v) {
        size_t k = (size_t)(i >> SH);
        if (k >= pages.size()) pages.resize(k + 1, nullptr);
        if (!pages[k]) pages[k] = (int*)calloc(PS, sizeof(int));
        pages[k][i & (PS - 1)] = v;
    }
};
static PagedArray hLenA, hBegA;
static set<int> ba[SM + 2];
static set<pair<int,int>> big;

static inline int bk(ll l) { return (int)(l < SM + 1 ? l : SM + 1); }
static inline ll hLen(ll s) { return hLenA.get(s); }
static inline void addRaw(ll s, ll len) {
    if (len <= 0) return;
    hLenA.set(s, (int)len); hBegA.set(s + len, (int)s + 1);
    ba[bk(len)].insert((int)s);
    if (len > SM) big.insert({(int)len, (int)s});
}
static inline void delRaw(ll s, ll len) {
    hLenA.set(s, 0); hBegA.set(s + len, 0);
    ba[bk(len)].erase((int)s);
    if (len > SM) big.erase({(int)len, (int)s});
}
static inline ll holeEndingAt(ll e) { int b = hBegA.get(e); return b ? b - 1 : -1; }
static void addFree(ll s, ll len) {           // освободить участок со склейкой соседей
    ll ps = holeEndingAt(s);
    if (ps >= 0) { ll pl = s - ps; delRaw(ps, pl); s = ps; len += pl; }
    ll nl = hLen(s + len);
    if (nl) { delRaw(s + len, nl); len += nl; }
    addRaw(s, len);
}

static bool leftFit(ll L, ll& hs, ll& hl) {   // самая левая подходящая дыра
    ll b = LLONG_MAX; int bl = -1;
    for (int k = (int)L; k <= SM + 1; k++)
        if (!ba[k].empty() && *ba[k].begin() < b) { b = *ba[k].begin(); bl = k; }
    if (bl < 0) return false;
    hs = b; hl = hLen(b); return true;
}
static bool rightFit(ll L, ll& hs, ll& hl) {  // самая правая подходящая дыра
    ll b = -1; int bl = -1;
    for (int k = (int)L; k <= SM + 1; k++)
        if (!ba[k].empty() && *ba[k].rbegin() > b) { b = *ba[k].rbegin(); bl = k; }
    if (bl < 0) return false;
    hs = b; hl = hLen(b); return true;
}

// Первая непустая корзина размеров [lo, hi] (по возрастанию), минуя дыру ban; из равных по
// размеру берём левую. Дыры длиннее 40 перебираем по возрастанию длины.
static bool firstBucket(ll lo, ll hi, ll ban, ll& hs, ll& hl) {
    for (ll k = lo; k <= hi && k <= SM; k++) {
        if (ba[k].empty()) continue;
        auto it = ba[k].begin();
        if (*it == ban && ++it == ba[k].end()) continue;
        hs = *it; hl = k; return true;
    }
    if (hi > SM) {
        for (auto it = big.begin(); it != big.end(); ++it) {
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
static bool goodFit(ll L, ll& hs, ll& hl) {
    ll ban = -1;
    if (!ba[(int)L].empty()) {                            // есть точное совпадение
        ll es = *ba[(int)L].begin();
        if (trenchLen > GUARD) { hs = es; hl = L; return true; }
        ll rs, rl;
        if (!rightFit(L, rs, rl) || es != rs) { hs = es; hl = L; return true; }
        ban = es;                                         // точная дыра самая правая - поищем другую
    }
    if (firstBucket(L + KMIN, LLONG_MAX, ban, hs, hl)) return true;
    if (firstBucket(L, L + KMIN - 1, ban, hs, hl)) return true;
    if (ban >= 0) { hs = ban; hl = L; return true; }
    return false;
}

// Наименьший период p (PRIOR_MIN_PERIOD..PRIOR_MAX_PERIOD), которому подчиняются последние PRIOR_REPS*p длин,
// при условии что в периоде есть разные длины.
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
// Ответ интерактору: числа через пробел и перевод строки, сразу со сбросом буфера
// (протокол пошаговый: следующую команду интерактор пришлёт только после ответа).
static inline char* putNum(char* p, ll v) {
    char t[24]; int n = 0;
    do { t[n++] = (char)('0' + v % 10); v /= 10; } while (v);
    while (n) *p++ = t[--n];
    return p;
}
static inline void answer2(ll a, ll b) {
    char buf[48], *p = putNum(buf, a); *p++ = ' '; p = putNum(p, b); *p++ = '\n';
    fwrite(buf, 1, (size_t)(p - buf), stdout); fflush(stdout);
}
static inline void answer1(ll a) {
    char buf[24], *p = putNum(buf, a); *p++ = '\n';
    fwrite(buf, 1, (size_t)(p - buf), stdout); fflush(stdout);
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
    vector<int> runLens;                       // длины укладок после последнего удаления (для режима 3)
    vector<signed char> predDoom;              // предсказание режима 3: -1 нет, 0 останется, 1 удалят
    predDoom.reserve(200005); predDoom.push_back(-1);
    ll priorHit = 0, priorMiss = 0; bool priorBad = false;
    vector<int> remIdx;                        // позиции удалённых труб текущей пачки
    bool havePattern = false, expTrusted = false;
    int patPeriod = 2, patMask = 0, patRunLen = 0, patRunLenPrev = 0, expLen = 0;   // шаблон: удаляются позиции i, у которых (patMask >> (i % patPeriod)) & 1
    int remBurst = 0;                          // удалений подряд после последней укладки
    // Режим 4: период по всем укладкам; судьба трубы известна через L4_H укладок.
    vector<int> allLens; allLens.reserve(200005);
    vector<signed char> lpPer, lpRes, lpPred;
    lpPer.reserve(200005); lpRes.reserve(200005); lpPred.reserve(200005);
    lpPer.push_back(0); lpRes.push_back(0); lpPred.push_back(-1);
    static ll jDoom[PRIOR_MAX_PERIOD + 1][PRIOR_MAX_PERIOD], jKeep[PRIOR_MAX_PERIOD + 1][PRIOR_MAX_PERIOD];
    int lastPer = 0, perStable = 0; ll l4Hit = 0, l4Miss = 0; bool l4Bad = false;

    for (;;) {
        int cmd = readInt();
        if (cmd == 0 || cmd == INT_MIN) break;

        if (cmd > 0) {
            ll L = cmd;

            bool burst = curCount >= MIN_RUN && remBurst >= 5 && remBurst * 8 >= curCount;
            remBurst = 0;
            if ((int)L != curSize || burst) {  // пачка сменилась (или закончилась всплеском удалений): анализируем предыдущую
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
            bool phase  = havePattern && patRunLen >= MIN_RUN;          // режим 2
            bool doomed = phase && (((patMask >> (idx % patPeriod)) & 1) != 0);
            bool prior = false;                // режим 3: решение по предположению
            bool learned = false;              // режим 4: решение по выученной судьбе позиций периода
            runLens.push_back((int)L);

            int lper = 0, lres = 0;
            ll l4dig = 0;                      // сколько метров обречённых труб придёт за период начиная с текущей
            allLens.push_back((int)L);
            {   // судим трубу, уложенную L4_H укладок назад: удалили её за это время или нет
                ll id = (ll)pipes.size() - L4_H;
                if (id >= 1 && lpPer[id] > 0) {
                    bool died = pipes[id].first < 0;
                    (died ? jDoom : jKeep)[lpPer[id]][lpRes[id]]++;
                    if (lpPred[id] >= 0) { if ((lpPred[id] == 1) == died) l4Hit++; else l4Miss++; }
                    if (l4Hit + l4Miss >= 20 && l4Miss * 5 > l4Hit + l4Miss) l4Bad = true;   // >20% ошибок - выключаемся
                }
            }
            lper = detectPeriod(allLens);
            perStable = (lper > 0 && lper == lastPer) ? perStable + 1 : 0;
            lastPer = lper;
            if (lper > 0) {
                lres = (int)((allLens.size() - 1) % lper);
                if (!phase && !l4Bad && perStable >= 2 * lper) {
                    // верим, если каждая позиция периода уже оценена >= 2 раз и удаляется либо почти
                    // всегда, либо почти никогда, и обречённые позиции вообще есть
                    bool pure = true, any = false;
                    for (int r = 0; r < lper && pure; r++) {
                        ll n = jDoom[lper][r] + jKeep[lper][r];
                        if (n < 2) { pure = false; break; }
                        if (jDoom[lper][r] * 100 > n * (100 - L4_PURE) && jDoom[lper][r] * 100 < n * L4_PURE) pure = false;
                        if (jDoom[lper][r] * 100 >= n * L4_PURE) any = true;
                    }
                    if (pure && any) {
                        learned = true; phase = true;
                        ll n = jDoom[lper][lres] + jKeep[lper][lres];
                        doomed = jDoom[lper][lres] * 2 >= n;
                        for (int j = 0; j < lper; j++) {
                            int r = (lres + j) % lper;
                            ll nn = jDoom[lper][r] + jKeep[lper][r];
                            if (jDoom[lper][r] * 2 >= nn) l4dig += allLens[allLens.size() - 1 + j - lper];
                        }
                    }
                }
            }
            if (!phase && !priorBad) {
                int per = detectPeriod(runLens);
                if (per > 0) {
                    int res = (int)(runLens.size() - 1) % per;
                    phase = true; prior = true;
                    doomed = (res % 2 == 0) && !(per % 2 == 1 && res == per - 1);
                }
            }

            ll hs = -1, hl = -1, D, dig = 0;
            bool got;
            if (learned && !doomed) {          // режим 4, остающаяся: точная дыра, иначе хвост траншеи
                got = !ba[(int)L].empty();
                if (got) { hs = *ba[(int)L].begin(); hl = L; }
            }
            else if (learned) got = goodFit(L, hs, hl);   // режим 4, обречённая: переиспользует дыры
            else if (!phase) got = goodFit(L, hs, hl);
            else got = doomed ? rightFit(L, hs, hl) : leftFit(L, hs, hl);
            if (got) {
                delRaw(hs, hl);
                if (!doomed || learned) { D = hs;           if (hl > L) addRaw(hs + L, hl - L); }
                else                    { D = hs + hl - L;  if (hl > L) addRaw(hs, hl - L); }
            } else {
                ll ts = trench, tl = 0;
                ll e = holeEndingAt(trench);
                if (e >= 0) { ts = e; tl = trench - ts; delRaw(ts, tl); }
                ll need = L - tl;
                if (learned && doomed && l4dig - tl > need) need = l4dig - tl;   // место под все обречённые трубы периода
                if (phase && !prior && !learned) {   // режим 2: копаем сразу на остаток пачки
                    ll rest = (ll)(expLen - idx) * L - tl;
                    ll cap = CAP_MULT * (ll)idx * L + 40;
                    if (!expTrusted && rest > cap) rest = cap;
                    if (rest > need) need = rest;
                }
                dig = need > 0 ? (need + 19) / 20 : 0;   // хвост может быть уже длиннее нужного (режим 4)
                trench += 20 * dig;
                ll nl = tl + 20 * dig;
                if (!doomed) { D = ts;             if (nl > L) addRaw(ts + L, nl - L); }
                else         { D = ts + nl - L;    if (nl > L) addRaw(ts, nl - L); }
            }

            pipes.push_back({D, L}); runIdx.push_back(idx); runId.push_back(curRun);
            predDoom.push_back(prior ? (doomed ? 1 : 0) : -1);
            lpPer.push_back((signed char)lper); lpRes.push_back((signed char)lres);
            lpPred.push_back(learned ? (doomed ? 1 : 0) : -1);
            answer2(dig, D);
        } else {
            ll k = -cmd;
            answer1(pipes[k].first);
            addFree(pipes[k].first, pipes[k].second);
            if (runId[k] == curRun) remIdx.push_back(runIdx[k]);
            remBurst++;
            if (predDoom[k] >= 0) {                   // сверяем предположение режима 3 с реальностью
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
