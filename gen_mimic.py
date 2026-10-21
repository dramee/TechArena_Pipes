# Генератор «мимиков» открытых тестов: та же структура и ловушка, что в открытом тесте N,
# но другие длины, размеры и случайные числа. Для каждого открытого теста - по 2 мимика:
# mNN_a.in, mNN_b.in (вес как у оригинала). Запуск: python gen_mimic.py -> папка mimic/
import random, math, os

class S:
    def __init__(s, seed):
        s.r = random.Random(seed); s.c = []; s.pid = 0; s.alive = []; s.pos = {}
    def lay(s, L):
        L = max(1, min(40, int(L))); s.pid += 1; s.c.append(L)
        s.pos[s.pid] = len(s.alive); s.alive.append(s.pid); return s.pid
    def rem(s, p):
        if p not in s.pos: return False
        i = s.pos.pop(p); last = s.alive.pop()
        if last != p: s.alive[i] = last; s.pos[last] = i
        s.c.append(-p); return True
    def by_rank(s, k):   # удалить k-ю по свежести живую трубу (0 - самая свежая)
        if not s.alive: return
        k = min(k, len(s.alive) - 1)
        if k < 400:
            top = sorted(s.alive[-min(len(s.alive), 4 * (k + 1) + 64):])
            if len(top) > k: s.rem(top[-1 - k]); return
        srt = sorted(s.alive); s.rem(srt[-1 - k])

# --- крошечные ручные (открытые 1-4, 6): чередование, удаление по чётности, потом трубы длиннее
def tiny(s, v):
    r = s.r
    if v == 1:     # как тест 1: 3 трубы, удаление, ещё одна, удаление
        a = [s.lay(r.randint(1, 20)) for _ in range(3)]; s.rem(a[1]); s.lay(r.randint(1, 10)); s.rem(a[0])
    elif v == 2:   # как тест 2: x 1 y 1 ..., удалить все единички, потом трубы растущей длины
        k = r.randint(3, 6); ids = []
        for i in range(k): ids.append(s.lay(1)); ids.append(s.lay(i + 2))
        for i in range(0, len(ids), 2): s.rem(ids[i + 1]) if r.random() < 0.0 else s.rem(ids[i])
        for i in range(k): s.lay(r.randint(1, k + 1))
    elif v == 3:   # как тест 3: несколько коротких, удаление нечётных, вставка
        ids = [s.lay(r.randint(2, 5)) for _ in range(r.randint(5, 8))]
        for p in ids[::2]: s.rem(p)
        s.lay(r.randint(1, 3)); s.rem(ids[-1]); s.lay(r.randint(4, 7))
    elif v == 4:   # как тест 4: две близкие длины вперемешку, удалить почти всё, затем длиннее
        a, b = r.randint(3, 8), 0; b = a + 1
        ids = [s.lay(r.choice([a, b])) for _ in range(r.randint(8, 12))]
        for p in r.sample(ids, len(ids) - 3): s.rem(p)
        for d in range(2, 5): s.lay(b + d - 1)
    else:          # как тест 6: мелкие + одна длинная, удаление нечётных
        ids = [s.lay(r.randint(1, 3)) for _ in range(5)] + [s.lay(r.randint(10, 15))]
        for p in ids[0:5:2]: s.rem(p)
        s.lay(1); s.lay(r.randint(2, 4))

# --- открытый 5: лестница k..k+m по n штук, удаляют каждую вторую (начиная с первой)
def ladder(s, start, end, n):
    for L in range(start, end + 1):
        ids = [s.lay(L) for _ in range(n)]
        for i, p in enumerate(ids):
            if i % 2 == 0: s.rem(p)
# --- открытый 9: удваивающиеся длины по n штук, удаляют половину
def doubling(s, start, n):
    L = start
    while L <= 40:
        ids = [s.lay(L) for _ in range(n)]
        for i, p in enumerate(ids):
            if i % 2 == 0: s.rem(p)
        L *= 2
# --- открытый 7: k 1 (k-1) 1 ... потом удаляют пары и кладут k, k+1
def t7like(s, top):
    ids = []
    for k in range(2, top + 1):
        ids.append(s.lay(k)); ids.append(s.lay(1)); ids.append(s.lay(k - 1)); ids.append(s.lay(1))
    j = 0
    for k in range(1, top):
        if j + 2 >= len(ids): break
        s.rem(ids[j]); s.rem(ids[j + 2]); j += 4
        s.lay(k); s.lay(k + 1)
# --- открытый 10: периодическая последовательность мелких длин; затем разом удаляют позиции периода
#     по маске (в оригинале - чётные, кроме последней), потом в каждом периоде удаляют ещё одну
#     оставшуюся трубу длины x и кладут x+1 и 1
def t10like(s, per_len, reps, odd=False):
    r = s.r; per = [r.randint(1, 5) for _ in range(per_len)]
    ids = []
    for _ in range(reps):
        for x in per: ids.append(s.lay(x))
    if odd: dm = [j % 2 == 1 for j in range(per_len)]
    else: dm = [j % 2 == 0 and not (per_len % 2 == 1 and j == per_len - 1) for j in range(per_len)]
    for t in range(reps):
        for j in range(per_len):
            if dm[j]: s.rem(ids[t * per_len + j])
    keep = [j for j in range(per_len) if not dm[j]]; j = r.choice(keep); x = per[j]
    for t in range(reps):
        s.rem(ids[t * per_len + j]); s.lay(min(40, x + 1)); s.lay(1)
# --- открытые 11 и 15: короткие 1..7 (20% удалений свежих), затем длинные 8..hi (50% удалений,
#     ранги с тяжёлым хвостом - удаляют и старые короткие)
def twophase(s, n, hi, split=0.5):
    r = s.r; n1 = int(n * split); t2 = 0
    for i in range(n):
        if i < n1:
            if s.alive and r.random() < 0.2:
                k = int(r.expovariate(1 / 22)) if r.random() < 0.95 else r.randint(0, 600)
                s.by_rank(k)
            else: s.lay(r.randint(1, 7))
        else:
            t2 += 1
            if r.random() < 0.5:
                if r.random() < 0.6: k = int(r.expovariate(1 / 45))
                else: k = r.randint(0, max(1, int(t2 * 0.35)))
                s.by_rank(k)
            else: s.lay(r.randint(8, hi) if r.random() > 0.01 else r.randint(1, 7))
# --- открытые 8, 12, 13: смесь «фрагментов»
def fragments(s, n, big20=False):
    r = s.r
    while len(s.c) < n:
        f = r.random()
        if f < 0.3:      # рост с заменой: кладём s, удаляем предыдущую трубу фрагмента
            a = r.randint(5, 12); b = min(40, a + r.randint(3, 10)); prev = None
            for L in range(a, b + 1):
                for _ in range(r.randint(8, 25)):
                    p = s.lay(L)
                    if prev: s.rem(prev)
                    prev = p
        elif f < 0.55:   # чередование двух длин (5/10, 5/13, ...) с удалением недавних вперемешку
            x, y = r.choice([(5, 10), (5, 13), (5, 15), (10, 20), (6, 12)]); buf = []
            for i in range(r.randint(40, 200)):
                buf.append(s.lay(x if i % 2 == 0 else y))
                if len(buf) > 4: s.rem(buf.pop(r.randrange(len(buf) - 2)))
        elif f < 0.7:    # серия одинаковых (20 или 10) с заменами
            L = 20 if (big20 or r.random() < 0.5) else 10; ids = []
            for i in range(r.randint(10, 40)):
                ids.append(s.lay(L))
                if r.random() < 0.4 and len(ids) > 3: s.rem(ids.pop(r.randrange(len(ids) - 1)))
        else:            # шум: случайные длины 5..40, изредка удаляют старые
            for i in range(r.randint(10, 60)):
                if s.alive and r.random() < 0.25: s.rem(r.choice(s.alive))
                else: s.lay(r.choice([5, 10, 15, 20, r.randint(5, 40)]))
# --- открытый 14: мелкие 1..3 (тройки через одну), длинная серия единичек без удалений,
#     единички с редкими длинными, в конце трубы 10..24 с удалениями
def t14like(s, n):
    r = s.r; a = int(n * 0.1); b = int(n * 0.6); c = int(n * 0.9)
    for i in range(n):
        if i < a: s.lay(3 if i % 2 == 0 else r.choice([1, 1, 2]))
        elif i < b: s.lay(1)
        elif i < c:
            if r.random() < 0.05 and s.alive: s.by_rank(r.choice([1, 2, 2, 7, r.randint(0, 3000)]))
            else: s.lay(1 if r.random() < 0.9 else r.randint(2, 39))
        else:
            if r.random() < 0.15 and s.alive: s.by_rank(int(r.expovariate(1 / 60)))
            else: s.lay(1 if r.random() < 0.45 else r.randint(10, 24))

PLAN = {   # открытый тест: (вес, [генератор мимика a, генератор мимика b])
    1: (0.1, [lambda s: tiny(s, 1), lambda s: tiny(s, 1)]),
    2: (1, [lambda s: tiny(s, 2), lambda s: tiny(s, 2)]),
    3: (1, [lambda s: tiny(s, 3), lambda s: tiny(s, 3)]),
    4: (1, [lambda s: tiny(s, 4), lambda s: tiny(s, 4)]),
    5: (1, [lambda s: ladder(s, 1, 40, 150), lambda s: ladder(s, 3, 38, 101)]),
    6: (1, [lambda s: tiny(s, 6), lambda s: tiny(s, 6)]),
    7: (1, [lambda s: t7like(s, 35), lambda s: t7like(s, 40)]),
    8: (1, [lambda s: fragments(s, 6500), lambda s: fragments(s, 7000)]),
    9: (1, [lambda s: doubling(s, 1, 800), lambda s: doubling(s, 3, 1200)]),
    10: (1, [lambda s: t10like(s, 7, 770), lambda s: t10like(s, 9, 600, odd=True)]),
    11: (1, [lambda s: twophase(s, 10000, 32), lambda s: twophase(s, 10000, 36, 0.45)]),
    12: (4, [lambda s: fragments(s, 1500, True), lambda s: fragments(s, 2500, True)]),
    13: (4, [lambda s: fragments(s, 10000), lambda s: fragments(s, 12000)]),
    14: (4, [lambda s: t14like(s, 25000), lambda s: t14like(s, 30000)]),
    15: (10, [lambda s: twophase(s, 200000, 40), lambda s: twophase(s, 180000, 40, 0.55)]),
}
if __name__ == '__main__':
    os.makedirs('mimic', exist_ok=True)
    for t, (w, gens) in PLAN.items():
        for j, g in enumerate(gens):
            s = S(7000 + 10 * t + j); g(s); c = s.c[:199999]
            fn = 'mimic/m%02d_%s.in' % (t, 'ab'[j])
            with open(fn, 'w', newline='') as f:
                f.write('%s %d\r\n' % (w, len(c))); f.write(''.join('%d\r\n' % x for x in c)); f.write('0\r\n')
            print(fn, len(c))
