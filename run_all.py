#!/usr/bin/env python3
"""Run the solution on every test through the interactor and save the total KPI.

Usage: python run_all.py [solution] [-d DIR ...] [--timeout SEC] [-o kpi_results]
Defaults: compiles solution.cpp with g++ -O2 (pass a binary to skip), test dirs: tests myTests.
Writes <out>.csv (per test) and <out>.json (summary with total KPI).
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = ".exe" if os.name == "nt" else ""


def test_key(path):
    name = os.path.basename(path)
    stem = name.split(".")[0]
    return (os.path.dirname(path), int(stem) if stem.isdigit() else 10**9, name)


def read_weight(test_file):
    with open(test_file) as f:
        return float(f.readline().split()[0])


def run_one(solution, interactor, test_file, timeout):
    """Runs in a scratch cwd because the interactor always writes ./score.txt."""
    with tempfile.TemporaryDirectory() as tmp:
        answer = os.path.join(tmp, "answer.txt")
        open(answer, "w").close()
        start = time.time()
        status = "ok"
        inter = subprocess.Popen(
            [interactor, test_file, "score.txt", answer, "report.txt"],
            cwd=tmp, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL)
        sol = subprocess.Popen(solution, cwd=tmp, stdin=inter.stdout, stdout=inter.stdin,
                               stderr=subprocess.DEVNULL)
        inter.stdout.close()  # parent must not hold pipe ends, or EOF never arrives
        try:
            sol.wait(timeout=timeout)
            inter.stdin.close()
            inter.wait(timeout=10)
        except subprocess.TimeoutExpired:
            status = "timeout"
            sol.kill()
            inter.kill()
        finally:
            sol.wait()
            inter.wait()
            if not inter.stdin.closed:
                inter.stdin.close()
        elapsed = time.time() - start
        score = 0.0
        try:
            with open(os.path.join(tmp, "score.txt")) as f:
                score = float(f.readline())
        except (OSError, ValueError):
            if status == "ok":
                status = "no_score"  # interactor rejected the output / crash
        return score, elapsed, status


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("solution", nargs="?", default=None,
                    help="prebuilt solution binary; if omitted, solution.cpp is compiled")
    ap.add_argument("--interactor", default=os.path.join(HERE, "interactor" + EXE))
    ap.add_argument("-d", dest="dirs", nargs="+", default=[os.path.join(HERE, "tests"), os.path.join(HERE, "myTests")])
    ap.add_argument("--timeout", type=float, default=60.0)
    ap.add_argument("-o", dest="out", default=os.path.join(HERE, "kpi_results"))
    args = ap.parse_args()

    if args.solution is None:
        build = tempfile.mkdtemp(prefix="sol_")
        solution = os.path.join(build, "solution" + EXE)
        src = os.path.join(build, "solution.cpp")  # copy: g++ may choke on non-ASCII paths
        shutil.copyfile(os.path.join(HERE, "solution.cpp"), src)
        print("Compiling solution.cpp ...")
        built = subprocess.run(["g++", "-O2", "-std=c++17", "-o", solution, src],
                               capture_output=True, text=True).returncode == 0
        if not built:
            prebuilt = os.path.join(HERE, "solution")
            if not os.path.exists(prebuilt):
                sys.exit("g++ failed and no prebuilt ./solution found; pass a binary explicitly")
            print("WARNING: g++ failed, using prebuilt ./solution (may be stale)")
            shutil.copyfile(prebuilt, solution)  # extensionless binary can't be spawned on Windows
    else:
        solution = os.path.abspath(args.solution)
    interactor = os.path.abspath(args.interactor)
    for p in (solution, interactor):
        if not os.path.exists(p):
            sys.exit(f"Not found: {p}")

    tests = []
    for d in args.dirs:
        if not os.path.isdir(d):
            sys.exit(f"Not a directory: {d}")
        tests += sorted((os.path.join(os.path.abspath(d), f) for f in os.listdir(d) if f.endswith(".in")),
                        key=test_key)

    rows = []
    total = total_weight = 0.0
    for t in tests:
        weight = read_weight(t)
        score, elapsed, status = run_one(solution, interactor, t, args.timeout)
        total += score
        total_weight += weight
        label = os.path.join(os.path.basename(os.path.dirname(t)), os.path.basename(t))
        rows.append([label, weight, f"{score:.9f}", f"{score / weight:.6f}", f"{elapsed:.2f}", status])
        print(f"{label:<14} w={weight:<5g} score={score:.6f} ratio={score / weight:.4f} "
              f"time={elapsed:.2f}s {status}", flush=True)

    kpi = {
        "solution": solution,
        "tests": len(rows),
        "failed": sum(r[5] != "ok" for r in rows),
        "total_score": total,
        "max_possible": total_weight,
        "kpi": total / total_weight if total_weight else 0.0,
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(args.out + ".csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["test", "weight", "score", "ratio", "time_s", "status"])
        w.writerows(rows)
    with open(args.out + ".json", "w") as f:
        json.dump(kpi, f, indent=2)
    print(f"\nTotal score: {total:.6f} / {total_weight:g}  (KPI = {kpi['kpi']:.6f}, failed: {kpi['failed']})")
    print(f"Saved: {args.out}.csv, {args.out}.json")


if __name__ == "__main__":
    main()
