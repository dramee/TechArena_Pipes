#!/usr/bin/env python3
"""Run the solution on every test through the interactor and save the total KPI.

Usage: python run_all.py [solution] [-d DIR ...] [--timeout SEC] [-o kpi_results]
Defaults: compiles solution.cpp with g++ -O2 (pass a binary to skip), test dirs: tests myTests.
Writes <out>.csv (per test, incl. solution CPU time) and <out>.json (summary with total KPI).
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


def cpu_seconds_windows(proc):
    """User+kernel CPU time of a finished process (Windows; handle is still open until Popen is freed)."""
    import ctypes
    from ctypes import wintypes
    ft = [wintypes.FILETIME() for _ in range(4)]
    ok = ctypes.windll.kernel32.GetProcessTimes(wintypes.HANDLE(int(proc._handle)), *[ctypes.byref(x) for x in ft])
    if not ok:
        return None
    val = lambda f: ((f.dwHighDateTime << 32) | f.dwLowDateTime) / 1e7
    return val(ft[2]) + val(ft[3])


def wait_with_cpu(proc, timeout):
    """Waits for proc; returns its CPU time in seconds (None if unknown). Raises TimeoutExpired."""
    if os.name == "nt":
        proc.wait(timeout=timeout)
        return cpu_seconds_windows(proc)
    deadline = time.time() + timeout
    while True:   # os.wait4 reaps the child and reports its own rusage
        pid, status, ru = os.wait4(proc.pid, os.WNOHANG)
        if pid:
            proc.returncode = os.waitstatus_to_exitcode(status)
            return ru.ru_utime + ru.ru_stime
        if time.time() > deadline:
            raise subprocess.TimeoutExpired(proc.args, timeout)
        time.sleep(0.005)


def run_one(solution, interactor, test_file, timeout):
    """Runs in a scratch cwd because the interactor always writes ./score.txt."""
    with tempfile.TemporaryDirectory() as tmp:
        answer = os.path.join(tmp, "answer.txt")
        open(answer, "w").close()
        start = time.time()
        status = "ok"
        cpu = None
        inter = subprocess.Popen(
            [interactor, test_file, "score.txt", answer, "report.txt"],
            cwd=tmp, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL)
        sol = subprocess.Popen(solution, cwd=tmp, stdin=inter.stdout, stdout=inter.stdin,
                               stderr=subprocess.DEVNULL)
        inter.stdout.close()  # parent must not hold pipe ends, or EOF never arrives
        try:
            cpu = wait_with_cpu(sol, timeout)
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
        return score, elapsed, cpu, status


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
    cpu_total = cpu_max = 0.0
    for t in tests:
        weight = read_weight(t)
        score, elapsed, cpu, status = run_one(solution, interactor, t, args.timeout)
        cpu_total += cpu or 0.0
        cpu_max = max(cpu_max, cpu or 0.0)
        total += score
        total_weight += weight
        label = os.path.join(os.path.basename(os.path.dirname(t)), os.path.basename(t))
        cpu_s = "?" if cpu is None else f"{cpu:.3f}"
        rows.append([label, weight, f"{score:.9f}", f"{score / weight:.6f}", f"{elapsed:.2f}", cpu_s, status])
        print(f"{label:<14} w={weight:<5g} score={score:.6f} ratio={score / weight:.4f} "
              f"time={elapsed:.2f}s cpu={cpu_s}s {status}", flush=True)

    kpi = {
        "solution": solution,
        "tests": len(rows),
        "failed": sum(r[6] != "ok" for r in rows),
        "solution_cpu_total_s": round(cpu_total, 3),
        "solution_cpu_max_s": round(cpu_max, 3),
        "total_score": total,
        "max_possible": total_weight,
        "kpi": total / total_weight if total_weight else 0.0,
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(args.out + ".csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["test", "weight", "score", "ratio", "time_s", "cpu_s", "status"])
        w.writerows(rows)
    with open(args.out + ".json", "w") as f:
        json.dump(kpi, f, indent=2)
    print(f"\nTotal score: {total:.6f} / {total_weight:g}  (KPI = {kpi['kpi']:.6f}, failed: {kpi['failed']})")
    print(f"Solution CPU: total {cpu_total:.2f}s, max per test {cpu_max:.3f}s")
    print(f"Saved: {args.out}.csv, {args.out}.json")


if __name__ == "__main__":
    main()
