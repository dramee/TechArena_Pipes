#!/usr/bin/env python3
"""Fuzz the solution through the real interactor: any rejected answer, crash or timeout is reported.

Usage: python fuzz.py [N] [SEED] [--solution BIN]
Generates N random tests from several families (random, staircase runs, periodic with removal
delay/noise/drift, phase changes, place-then-remove) and runs each through ./interactor.
Without --solution, solution.cpp is compiled. Failing tests are kept in the printed directory.
"""
import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_all  # noqa: E402


class B:
    def __init__(self): self.ops = []; self.n = 0; self.alive = []
    def place(self, l): self.ops.append(l); self.n += 1; self.alive.append(self.n); return self.n
    def rem(self, i): self.ops.append(-i); self.alive.remove(i)
    def rem_random(self, rng):
        j = rng.randrange(len(self.alive)); i = self.alive[j]
        self.alive[j] = self.alive[-1]; self.alive.pop(); self.ops.append(-i)


def fam_random(rng):
    b = B(); n = rng.choice([rng.randint(1, 30), rng.randint(30, 600), rng.randint(600, 4000)])
    lo = rng.randint(1, 40); hi = rng.randint(lo, 40); p = rng.random() * 0.7
    for _ in range(n):
        if b.alive and rng.random() < p:
            if rng.random() < 0.3: b.rem(b.alive[-1])
            else: b.rem_random(rng)
        else: b.place(rng.randint(lo, hi))
    return b.ops


def fam_staircase(rng):
    b = B(); s = rng.randint(1, 10); m = rng.randint(20, 400); per = rng.randint(2, 6)
    mask = rng.sample(range(per), rng.randint(1, per - 1))
    for _ in range(rng.randint(2, 10)):
        ids = [b.place(min(40, s)) for _ in range(m if rng.random() < 0.7 else rng.randint(5, 2 * m))]
        dead = [ids[i] for i in range(len(ids)) if i % per in mask]
        if rng.random() < 0.3: rng.shuffle(dead)
        for d in dead: b.rem(d)
        s += rng.choice([0, 1, 1, 2, 5])
    return b.ops


def fam_periodic(rng):
    b = B(); per = rng.randint(2, 13); pat = [rng.randint(1, rng.choice([5, 12, 40])) for _ in range(per)]
    kill = set(rng.sample(range(per), rng.randint(0, per - 1))); delay = rng.randint(0, 3)
    noise = rng.choice([0, 0, 0.02, 0.2]); pend = []
    for _ in range(rng.randint(5, 600)):
        ids = [b.place(pat[i]) for i in range(per)]
        pend.append([ids[i] for i in range(per) if (i in kill) != (rng.random() < noise)])
        if len(pend) > delay:
            for d in pend.pop(0): b.rem(d)
        if rng.random() < 0.02: pat[rng.randrange(per)] = rng.randint(1, 40)   # pattern drift
        if rng.random() < 0.01: b.place(rng.randint(1, 40))                      # stray pipe
    return b.ops


def fam_phases(rng):
    b = B()
    for _ in range(rng.randint(2, 5)):
        lo = rng.randint(1, 40); hi = rng.randint(lo, 40); p = rng.random() * 0.8
        for _ in range(rng.randint(50, 1500)):
            if b.alive and rng.random() < p: b.rem_random(rng)
            else: b.place(rng.randint(lo, hi))
    return b.ops


def fam_place_then_kill(rng):
    b = B(); per = rng.randint(3, 9); pat = [rng.randint(1, 8) for _ in range(per)]
    ids = [b.place(pat[i % per]) for i in range(rng.randint(50, 3000))]
    for i in ids:
        if rng.random() < 0.5: b.rem(i)
        if rng.random() < 0.5: b.place(rng.randint(1, 40))
    return b.ops


FAMS = [fam_random, fam_staircase, fam_periodic, fam_phases, fam_place_then_kill]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("n", nargs="?", type=int, default=300)
    ap.add_argument("seed", nargs="?", type=int, default=1)
    ap.add_argument("--solution", default=None)
    ap.add_argument("--interactor", default=os.path.join(HERE, "interactor" + run_all.EXE))
    args = ap.parse_args()

    out = tempfile.mkdtemp(prefix="fuzz_")
    sol = args.solution
    if sol is None:
        sol = os.path.join(out, "solution" + run_all.EXE)
        src = os.path.join(out, "solution.cpp"); shutil.copyfile(os.path.join(HERE, "solution.cpp"), src)
        if subprocess.run(["g++", "-O2", "-std=c++17", "-o", sol, src]).returncode:
            sys.exit("g++ failed")
    sol, inter = os.path.abspath(sol), os.path.abspath(args.interactor)

    def job(k):
        rng = random.Random(args.seed * 100000 + k)
        fam = FAMS[k % len(FAMS)]
        ops = fam(rng) or [1]
        if not any(o > 0 for o in ops): ops.append(1)
        path = os.path.join(out, f"{fam.__name__}_{k}.in")
        with open(path, "w") as f:
            f.write(f"1 {len(ops)}\n" + "\n".join(map(str, ops)) + "\n0\n")
        score, _, _, status = run_all.run_one(sol, inter, path, 60)
        if status == "ok": os.remove(path)
        return path, status, len(ops)

    with ThreadPoolExecutor(6) as ex:
        res = list(ex.map(job, range(args.n)))
    bad = [r for r in res if r[1] != "ok"]
    print(f"{len(res)} tests ({sum(r[2] for r in res)} ops), {len(bad)} failed")
    for p, st, n in bad[:20]:
        print(" FAIL", st, p, n, "ops")
    if not bad: shutil.rmtree(out, ignore_errors=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
