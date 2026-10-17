// Gas pipeline construction. Онлайн-решение.
//
// KPI = сумма длин оставшихся труб / длина траншеи. Числитель задан сценарием,
// значит задача = минимизировать выкопанное. Это динамическое распределение памяти.
//
// Режим 1 (обычный): best-fit с приоритетом точного совпадения по размеру.
//   Порог KMIN: предпочитаем дыру, оставляющую остаток >= KMIN, чтобы не плодить
//   огрызки по 1-3 м. На коротких тестах лучше KMIN=4, на длинных KMIN=16,
//   поэтому после 500-й трубы порог переключается.
//
// Режим 2 (фазовый): включается, если обнаружена пачка из >= MIN_RUN труб одной
//   длины, в которой затем удалили каждую вторую. Тогда в следующих пачках:
//     - трубы, которые по выученной чётности будут удалены, кладём к концу траншеи;
//     - остальные - к началу;
//     - копаем сразу на всю оставшуюся пачку.
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
static const ll SWITCH_AT = 500;   // после скольких труб поднять порог остатка
static ll KMIN = 4;
static const ll GUARD = 100;   // правило правой дыры работает только на коротких траншеях
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
static bool goodFit(ll L, ll& hs, ll& hl) {
    ll rs, rl;
    bool anyRight = rightFit(L, rs, rl);                  // самая правая подходящая дыра
    if (!ba[(int)L].empty()) {                            // есть точное совпадение
        ll es = *ba[(int)L].begin();
        if (!anyRight || es != rs || trenchLen > GUARD) { hs = es; hl = L; return true; }
        if (pickExcept(L, es, hs, hl)) return true;       // точная дыра самая правая — ищем другую
        hs = es; hl = L; return true;                     // других нет, берём её
    }
    auto it = bySize.lower_bound({L + KMIN, LLONG_MIN});
    if (it == bySize.end()) it = bySize.lower_bound({L, LLONG_MIN});
    if (it == bySize.end()) return false;
    hl = it->first; hs = it->second; return true;
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
    vector<pair<ll,ll>> pipes;  pipes.reserve(200005);  pipes.push_back({-1,-1});
    vector<int> runIdx;         runIdx.reserve(200005); runIdx.push_back(-1);
    vector<int> runId;          runId.reserve(200005);  runId.push_back(-1);

    int curSize = -1, curRun = 0, curCount = 0;
    vector<int> remIdx;                        // позиции удалённых труб текущей пачки
    bool havePattern = false;
    int patParity = 0, patRunLen = 0;

    for (;;) {
        int cmd = readInt();
        if (cmd == 0 || cmd == INT_MIN) break;

        if (cmd > 0) {
            ll L = cmd;

            if ((int)L != curSize) {           // пачка сменилась: анализируем предыдущую
                if (curSize >= 0 && curCount >= MIN_RUN && !remIdx.empty()) {
                    bool same = true; int par = remIdx[0] & 1;
                    for (size_t i = 1; i < remIdx.size(); i++)
                        if ((remIdx[i] & 1) != par) { same = false; break; }
                    ll cnt = (ll)remIdx.size();
                    if (same && cnt * 3 >= curCount && cnt * 3 <= (ll)curCount * 2) {
                        havePattern = true; patParity = par; patRunLen = curCount;
                    } else havePattern = false;
                }
                curSize = (int)L; curRun++; curCount = 0; remIdx.clear();
            }

            if (++laid == SWITCH_AT) KMIN = 16;
            int idx = curCount++;
            bool phase  = havePattern && patRunLen >= MIN_RUN;
            bool doomed = phase && ((idx & 1) == patParity);

            ll hs = -1, hl = -1, D, dig = 0;
            bool got = !phase ? goodFit(L, hs, hl)
                              : (doomed ? rightFit(L, hs, hl) : leftFit(L, hs, hl));
            if (got) {
                delRaw(hs, hl);
                if (!doomed) { D = hs;            if (hl > L) addRaw(hs + L, hl - L); }
                else         { D = hs + hl - L;   if (hl > L) addRaw(hs, hl - L); }
            } else {
                ll ts = trench, tl = 0;
                auto e = fe.find(trench);
                if (e != fe.end()) { ts = e->second; tl = trench - ts; delRaw(ts, tl); }
                ll need = L - tl;
                if (phase) {                      // копаем сразу на остаток пачки
                    ll rest = (ll)(patRunLen - idx) * L - tl;
                    if (rest > need) need = rest;
                }
                dig = (need + 19) / 20;
                trench += 20 * dig;
                ll nl = tl + 20 * dig;
                if (!doomed) { D = ts;             if (nl > L) addRaw(ts + L, nl - L); }
                else         { D = ts + nl - L;    if (nl > L) addRaw(ts, nl - L); }
            }

            pipes.push_back({D, L}); runIdx.push_back(idx); runId.push_back(curRun);
            printf("%lld %lld\n", dig, D);
            fflush(stdout);
        } else {
            ll k = -cmd;
            printf("%lld\n", pipes[k].first);
            fflush(stdout);
            addFree(pipes[k].first, pipes[k].second);
            if (runId[k] == curRun) remIdx.push_back(runIdx[k]);
            pipes[k] = {-1, -1};
        }
    }
    return 0;
}
