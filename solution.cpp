#include <cstdio>
#include <algorithm>
#include <iterator>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

static const double EPS = 0.5;
static const long long SWITCH_AT = 100000;    // до этого числа труб включительно: среднее, после: медиана

static unordered_map<int,int> freeByStart;  // start -> end
static unordered_map<int,int> freeByEnd;    // end -> start
static set<int> buckets[42];                // дырки размера 1..40, корзина 41 = больше 40
static set<pair<int,int>> bySize;           // (размер, start) всех дырок

static long long lenCount[41];              // сколько труб каждой длины уже получено
static long long cnt = 0;                   // всего поступило труб
static long long sumLen = 0;                // сумма длин поступивших труб

static inline void addFree(int s, int e) {
    freeByStart[s] = e;
    freeByEnd[e] = s;
    buckets[min(e - s, 41)].insert(s);
    bySize.insert({e - s, s});
}

static inline int removeFree(int s) {
    int e = freeByStart[s];
    freeByStart.erase(s);
    freeByEnd.erase(e);
    buckets[min(e - s, 41)].erase(s);
    bySize.erase({e - s, s});
    return e;
}

// значение с порядковым номером k (с нуля) в отсортированной последовательности длин
static int kthLength(long long k) {
    long long acc = 0;
    for (int v = 1; v <= 40; v++) {
        acc += lenCount[v];
        if (k < acc) return v;
    }
    return 40;
}

static double median() {
    int a = kthLength((cnt - 1) / 2);
    int b = kthLength(cnt / 2);
    return (a + b) / 2.0;
}

static double mean() {
    return (double)sumLen / (double)cnt;
}

// текущая статистика: среднее до SWITCH_AT труб, затем медиана (cnt > 0)
static double stat() {
    return cnt <= SWITCH_AT ? mean() : median();
}

int main() {
    freeByStart.reserve(1 << 19);
    freeByEnd.reserve(1 << 19);

    vector<pair<int,int>> pipes(1);  // pid -> (start, length)
    int trench = 0;
    int cmd;

    while (scanf("%d", &cmd) == 1) {
        if (cmd == 0) break;

        if (cmd > 0) {
            int L = cmd;
            int s = -1;

            // труба "маленькая" относительно среднего/медианы по предыдущим трубам
            bool isSmall = (cnt > 0) && ((double)L < EPS * stat());

            // маленькую трубу пробуем положить в самую большую дырку
            if (isSmall && !bySize.empty()) {
                auto it = prev(bySize.end());
                if (it->first >= L) s = it->second;
            }

            // иначе (или если не влезло) обычный best-fit
            if (s < 0) {
                for (int b = L; b < 42; b++) {
                    if (!buckets[b].empty()) {
                        s = *buckets[b].begin();
                        break;
                    }
                }
            }

            if (s >= 0) {
                int e = removeFree(s);
                if (s + L < e) addFree(s + L, e);
                pipes.push_back({s, L});
                printf("0 %d\n", s);
            } else {
                int ts;
                auto it = freeByEnd.find(trench);
                if (it == freeByEnd.end()) {
                    ts = trench;
                } else {
                    ts = it->second;
                    removeFree(ts);
                }
                int need = L - (trench - ts);
                int n = (need + 19) / 20;
                trench += 20 * n;
                if (ts + L < trench) addFree(ts + L, trench);
                pipes.push_back({ts, L});
                printf("%d %d\n", n, ts);
            }
            fflush(stdout);

            lenCount[L]++;
            sumLen += L;
            cnt++;
        } else {
            int k = -cmd;
            int s = pipes[k].first, L = pipes[k].second;
            printf("%d\n", s);
            fflush(stdout);

            int ns = s, ne = s + L;
            auto li = freeByEnd.find(s);
            if (li != freeByEnd.end()) {
                int ls = li->second;
                removeFree(ls);
                ns = ls;
            }
            if (freeByStart.count(ne)) {
                ne = removeFree(ne);
            }
            addFree(ns, ne);
        }
    }
    return 0;
}