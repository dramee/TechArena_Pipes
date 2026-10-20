"""Edge-case tests for the pipes solution -> edgeTests/*.in

Each test targets a specific place where solution.cpp could break (invalid answer, crash,
time/memory) or where one of its heuristics could make things much worse than plain best-fit.
Format: "<weight> <n>", n operations (L > 0 places a pipe, -k removes pipe k), then 0.
"""
import os
import random

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "edgeTests")


class B:
    """Operation builder that keeps track of alive pipes so every removal is valid."""
    def __init__(self):
        self.ops, self.n, self.alive = [], 0, set()

    def place(self, length):
        assert 1 <= length <= 40
        self.ops.append(length); self.n += 1; self.alive.add(self.n)
        return self.n

    def rem(self, pid):
        assert pid in self.alive
        self.ops.append(-pid); self.alive.remove(pid)

    def rem_random(self, rng):
        self.rem(rng.choice(sorted(self.alive)))


TESTS = []   # (name, weight, description, builder)


def test(name, weight, desc):
    def deco(fn):
        TESTS.append((name, weight, desc, fn)); return fn
    return deco


# ---------- tiny / boundary ----------
@test("e01_single1", 1, "одна труба длины 1")
def _(b, rng): b.place(1)


@test("e02_single40", 1, "одна труба длины 40 (две копки сразу)")
def _(b, rng): b.place(40)


@test("e03_place_remove_place", 1, "положить, удалить, положить снова: траншея одной дырой")
def _(b, rng):
    a = b.place(7); b.rem(a); b.place(13)


@test("e04_bucket_41", 1, "дыры ровно 40, 41 и 42 м: граница корзины «больше 40»")
def _(b, rng):
    ids = [b.place(x) for x in [40, 1, 40, 2, 40, 1, 39, 1, 40] * 30]
    for i, pid in enumerate(ids):
        if i % 9 in (0, 1, 2, 3): b.rem(pid)     # 40+1 = 41, 40+2 = 42 holes after merges
    for x in [40, 40, 39, 1, 2, 40] * 40: b.place(x)


@test("e05_min_run_edge", 1, "пачки ровно по 19/20/21 труб (граница MIN_RUN) с удалением каждой второй")
def _(b, rng):
    size = 1
    for _ in range(12):
        for cnt in (19, 20, 21):
            ids = [b.place(size) for _ in range(cnt)]
            for i in range(0, cnt, 2): b.rem(ids[i])
            size = min(40, size + 1)


@test("e06_burst_edge", 1, "одна длина подряд, всплески ровно 4 и 5 удалений (граница разбиения пачки)")
def _(b, rng):
    for rep in range(200):
        ids = [b.place(6) for _ in range(40)]
        k = 4 if rep % 2 else 5
        for pid in ids[:2 * k:2]: b.rem(pid)


# ---------- growth / structure extremes ----------
@test("e07_all40_200k", 10, "200 000 труб по 40 м без удалений: траншея 8 млн м")
def _(b, rng):
    for _ in range(200000): b.place(40)


@test("e08_all1_200k", 10, "200 000 труб по 1 м без удалений")
def _(b, rng):
    for _ in range(200000): b.place(1)


@test("e09_lifo_stack", 10, "стек: удаляется всегда последняя уложенная труба")
def _(b, rng):
    stack = []
    for _ in range(150000):
        if stack and rng.random() < 0.48: b.rem(stack.pop())
        else: stack.append(b.place(rng.randint(1, 40)))


@test("e10_fifo_queue", 10, "очередь: удаляется всегда самая старая труба")
def _(b, rng):
    q = []; head = 0
    for _ in range(150000):
        if len(q) - head > 500 and rng.random() < 0.5: b.rem(q[head]); head += 1
        else: q.append(b.place(rng.randint(1, 40)))


@test("e11_remove_all_refill", 4, "волны: заполнить, удалить всё, заполнить заново")
def _(b, rng):
    for wave in range(8):
        ids = [b.place(rng.randint(1, 40)) for _ in range(2000)]
        if wave < 7:
            rng.shuffle(ids)
            for pid in ids: b.rem(pid)


@test("e12_place_remove_instantly", 4, "каждую трубу удаляют сразу после укладки")
def _(b, rng):
    for _ in range(20000):
        b.rem(b.place(rng.randint(1, 40)))
    b.place(5)


@test("e13_sparse_survivors", 4, "удалить все, кроме каждой сотой, потом заполнять случайными")
def _(b, rng):
    ids = [b.place(rng.randint(1, 40)) for _ in range(10000)]
    for i, pid in enumerate(ids):
        if i % 100: b.rem(pid)
    for _ in range(10000): b.place(rng.randint(1, 40))


@test("e14_only_old_removed", 10, "удаляют только трубы старше 10 000 операций")
def _(b, rng):
    from collections import deque
    young, old = deque(), []          # young: (time, pid) younger than 10 000 ops; old: removable pool
    for t in range(150000):
        while young and t - young[0][0] > 10000: old.append(young.popleft()[1])
        if old and rng.random() < 0.4:
            j = rng.randrange(len(old)); old[j], old[-1] = old[-1], old[j]
            b.rem(old.pop())
        else: young.append((t, b.place(rng.randint(1, 40))))


@test("e15_all_removed_end", 1, "в конце удаляют все трубы (KPI 0 неизбежен, проверка на падение)")
def _(b, rng):
    ids = [b.place(rng.randint(1, 40)) for _ in range(500)]
    for pid in ids: b.rem(pid)
    b.place(3)


# ---------- phase mode (2): staircase runs ----------
def staircase(b, sizes, counts, removed):
    for s, c in zip(sizes, counts):
        ids = [b.place(s) for _ in range(c)]
        for i, pid in enumerate(ids):
            if removed(i): b.rem(pid)


@test("e16_trusted_then_short", 1, "режим 2: три одинаковые пачки по 2000 (прогноз подтверждён), затем короткая пачка")
def _(b, rng):
    staircase(b, [1, 2, 3, 4], [2000, 2000, 2000, 30], lambda i: i % 2 == 0)


@test("e17_trusted_then_kept", 1, "режим 2: шаблон подтверждён, потом удаления прекращаются")
def _(b, rng):
    staircase(b, [1, 2, 3], [1000, 1000, 1000], lambda i: i % 2 == 0)
    for s in (4, 5, 6): [b.place(s) for _ in range(1000)]


@test("e18_trusted_then_other_mask", 1, "режим 2: шаблон «каждая вторая» сменяется на «каждая третья»")
def _(b, rng):
    staircase(b, [1, 2, 3, 4], [900] * 4, lambda i: i % 2 == 0)
    staircase(b, [5, 6, 7, 8], [900] * 4, lambda i: i % 3 == 0)


@test("e19_full_run_removed", 1, "режим 2: пачки, удалённые целиком, вперемешку с половинными")
def _(b, rng):
    for k in range(1, 21):
        ids = [b.place(k) for _ in range(200)]
        for i, pid in enumerate(ids):
            if k % 3 == 0 or i % 2 == 0: b.rem(pid)


@test("e20_staircase_down", 1, "лестница вниз: 40, 39, ... 1 с удалением каждой второй")
def _(b, rng):
    staircase(b, list(range(40, 0, -1)), [150] * 40, lambda i: i % 2 == 1)


# ---------- prior mode (3): guessed periodic fate ----------
@test("e21_prior_wrong", 1, "режим 3: период 7 без удалений, потом удаляют нечётные позиции (догадка неверна)")
def _(b, rng):
    pat = [3, 5, 2, 5, 1, 3, 1]
    ids = [b.place(pat[i % 7]) for i in range(5000)]
    for i, pid in enumerate(ids):
        if (i % 7) % 2 == 1: b.rem(pid)
    for _ in range(1500): b.place(4)


@test("e22_prior_none_removed", 1, "режим 3: периодичный поток, удалений не бывает вовсе")
def _(b, rng):
    pat = [6, 2, 9, 2, 6]
    for i in range(20000): b.place(pat[i % 5])


# ---------- learned periodic mode (4) ----------
def periodic(b, rng, pat, kill, blocks, delay=0, stray_every=0, noise=0.0):
    pend = []
    for r in range(blocks):
        if stray_every and r and r % stray_every == 0:
            b.place(rng.randint(1, 40))                 # stray pipe: shifts global position in the period
        ids = [b.place(pat[i]) for i in range(len(pat))]
        pend.append([ids[i] for i in range(len(pat)) if (i in kill) != (rng.random() < noise)])
        while len(pend) > delay:
            for pid in pend.pop(0): b.rem(pid)


@test("e23_l4_switch_same_period", 4, "режим 4: шаблон и набор удаляемых позиций меняются, период тот же (5)")
def _(b, rng):
    periodic(b, rng, [8, 3, 5, 3, 8], {0, 2}, 2000, delay=1)
    periodic(b, rng, [4, 9, 4, 9, 4], {1, 3, 4}, 2000, delay=1)


@test("e24_l4_stray_shift", 4, "режим 4: каждые 37 блоков лишняя труба сдвигает позиции периода")
def _(b, rng):
    periodic(b, rng, [7, 2, 6, 2, 7, 5], {1, 3}, 3000, delay=1, stray_every=37)


@test("e25_l4_stray_shift_kept", 4, "режим 4: как e24, но лишняя труба есть в каждом 3-м блоке")
def _(b, rng):
    periodic(b, rng, [7, 2, 6, 2, 7, 5], {1, 3}, 3000, delay=0, stray_every=3)


@test("e26_l4_lag_boundary", 4, "режим 4: трубы удаляют примерно через 64 укладки (граница суждения L4_H)")
def _(b, rng):
    periodic(b, rng, [5, 9, 5, 9, 5, 9, 5, 9], {1, 3, 5}, 2000, delay=8)


@test("e27_l4_noise30", 4, "режим 4: 30% шума в судьбе труб")
def _(b, rng):
    periodic(b, rng, [6, 6, 2, 9, 3], {0, 3}, 3000, delay=1, noise=0.3)


@test("e28_period2", 1, "период 2 (меньше минимального 3): чередование 5 и 6, удаляют шестёрки")
def _(b, rng):
    periodic(b, rng, [5, 6], {1}, 5000, delay=1)


@test("e29_period13", 1, "период 13 (больше максимального 12)")
def _(b, rng):
    periodic(b, rng, [rng.randint(1, 12) for _ in range(13)], {0, 4, 7, 11}, 1500, delay=1)


@test("e30_almost_periodic_200k", 10, "194 000 укладок почти периодичны (период 12, каждый 12-й блок сбит): цена поиска периода")
def _(b, rng):
    pat = [rng.randint(1, 40) for _ in range(12)]
    for r in range(194000 // 12):
        blk = list(pat)
        if r % 12 == 11: blk[-1] = blk[-1] % 40 + 1
        for x in blk: b.place(x)
    for pid in sorted(b.alive)[::3][:5000]: b.rem(pid)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "README.txt"), "w", encoding="utf-8") as readme:
        for idx, (name, weight, desc, fn) in enumerate(TESTS):
            rng = random.Random(1000 + idx)
            b = B()
            fn(b, rng)
            if not any(o > 0 for o in b.ops): b.place(1)
            with open(os.path.join(OUT_DIR, name + ".in"), "w") as f:
                f.write(f"{weight} {len(b.ops)}\n" + "\n".join(map(str, b.ops)) + "\n0\n")
            readme.write(f"{name}.in  w={weight:<3} {len(b.ops):>7} ops  {desc}\n")
    print(f"Generated {len(TESTS)} tests in {OUT_DIR}")


if __name__ == "__main__":
    main()
