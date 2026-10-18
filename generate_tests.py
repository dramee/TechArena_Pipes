"""Local generator of tests (regular 16..65, stress 66..85) for the pipes problem (no external API).

Format: "<weight> <n>", then n operations (L in 1..40 places a pipe, -k removes
the k-th placed pipe that is still alive), then a terminating 0.
"""
import os
import random

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "myTests")
START_NUM = 16
COUNT = 70          # 16.in .. 85.in
STRESS_FROM = 66    # 16..65 - regular tests, 66..85 - stress tests
MAX_LEN = 40


def gen(n, rng, lo, hi, remove_p, lifo_p):
    """n operations; lengths in [lo, hi]; remove with prob remove_p,
    taking the newest alive pipe with prob lifo_p, else a random one."""
    ops, alive, next_id = [], [], 1
    for _ in range(n):
        if alive and rng.random() < remove_p:
            i = len(alive) - 1 if rng.random() < lifo_p else rng.randrange(len(alive))
            alive[i], alive[-1] = alive[-1], alive[i]
            ops.append(-alive.pop())
        else:
            ops.append(rng.randint(lo, hi))
            alive.append(next_id)
            next_id += 1
    return ops


def periodic(n, rng, period):
    """Repeating block of lengths with every pipe of the block removed afterwards."""
    ops, next_id = [], 1
    while len(ops) < n:
        block = [rng.randint(1, MAX_LEN) for _ in range(period)]
        ids = list(range(next_id, next_id + period))
        next_id += period
        ops += block
        rng.shuffle(ids)
        ops += [-i for i in ids[: period // 2]]
        # leftover half stays alive
    return ops[:n]


def make(num):
    rng = random.Random(num)
    if num >= STRESS_FROM:
        n = rng.randint(100000, 200000)
    else:
        n = int(10 ** rng.uniform(0.7, 4.4))      # ~5 .. 25000 operations
    kind = (num - START_NUM) % 5
    if kind == 0:
        ops = gen(n, rng, 1, MAX_LEN, 0.45, 0.2)
    elif kind == 1:
        ops = gen(n, rng, 1, 5, 0.48, 0.5)
    elif kind == 2:
        ops = gen(n, rng, 20, MAX_LEN, 0.4, 0.0)
    elif kind == 3:
        ops = periodic(n, rng, rng.choice([5, 10, 50, 200]))
    else:
        ops = gen(n // 2, rng, 1, MAX_LEN, 0.0, 0.0)
        alive = list(range(1, len(ops) + 1))
        rng.shuffle(alive)
        while len(ops) < n:
            if alive and rng.random() < 0.7:
                ops.append(-alive.pop())
            else:
                ops.append(rng.randint(1, MAX_LEN))
    return ops


def weight(num, n):
    if num >= STRESS_FROM:
        return 10
    return 1 if n < 10000 else 4


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for num in range(START_NUM, START_NUM + COUNT):
        ops = make(num)
        with open(os.path.join(OUT_DIR, f"{num}.in"), "w") as f:
            f.write(f"{weight(num, len(ops))} {len(ops)}\n")
            f.write("\n".join(map(str, ops)))
            f.write("\n0\n")
    print(f"Generated {COUNT} tests: {START_NUM}.in .. {START_NUM + COUNT - 1}.in")


if __name__ == "__main__":
    main()
