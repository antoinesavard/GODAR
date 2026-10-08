#!/usr/bin/env python3
"""
GODAR regression tests.

Runs the test cases described in tests/cases.json with a compiled godar
executable, then checks the outputs:

    reference  outputs of a run against the stored references
               (tests/reference/<case>/)
    same       outputs of two runs of the same case, e.g. 1 vs 2 MPI ranks
    restart    first output of a restarted run against the output of the
               run it restarted from
    mass       ice mass (sum of r**2 h over the floes) conserved in a run

Each check compares the files with a named set of tolerances from
cases.json. The run times are printed and appended to <work>/timings.csv.

Only the Python standard library is needed. From the repository root:

    python tests/run_tests.py                      # all cases
    python tests/run_tests.py --case TEST          # one case
    python tests/run_tests.py --update-reference   # store new references
"""

import argparse
import csv
import datetime
import gzip
import json
import math
import re
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

TESTS = Path(__file__).resolve().parent
ROOT = TESTS.parent

# outputs with one line of n values per output step
DENSE = ["x", "y", "u", "v", "theta", "omega", "r", "h", "tfx", "tfy",
         "mom", "tsigxx", "tsigyy", "tsigxy", "tsigyx", "tp"]

# outputs with one line "idx j i [values]" per bond or contact
SPARSE = ["bond", "damage", "beam", "angle"]


# ----------------------------------------------------------------------
#   reading the outputs
# ----------------------------------------------------------------------

def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def read_output(path, name):
    """dense: list of output steps, each a list of n values
    sparse: dict (idx, j, i) -> list of values"""
    with open_text(path) as f:
        if name in DENSE:
            return [[float(v) for v in line.split()]
                    for line in f if line.strip()]
        out = {}
        for line in f:
            p = line.split()
            if p:
                out[(int(p[0]), int(p[1]), int(p[2]))] = \
                    [float(v) for v in p[3:]]
        return out


def store_reference(src, dst):
    """gzip src into dst, without time stamp so that git sees no change
    when the content is the same"""
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(src, "rb") as fin, open(dst, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fout:
            shutil.copyfileobj(fin, fout)


# ----------------------------------------------------------------------
#   comparisons
# ----------------------------------------------------------------------

def max_diff(pairs):
    """largest |a - b| over the pairs of values, inf if a NaN is found"""
    worst = 0.0
    for a, b in pairs:
        d = abs(a - b)
        if not d <= worst:
            worst = d if d == d else math.inf
    return worst


def compare(new, ref, name, tol):
    """compare two outputs of the same file. tol = [atol, rtol]: the
    largest difference must stay below atol + rtol * max|ref|"""
    atol, rtol = tol
    if name in DENSE:
        if len(new) != len(ref) or \
                any(len(a) != len(b) for a, b in zip(new, ref)):
            return False, "different number of outputs or particles"
        ref_values = [v for line in ref for v in line]
        pairs = zip((v for line in new for v in line), ref_values)
    else:
        missing, extra = ref.keys() - new.keys(), new.keys() - ref.keys()
        if missing or extra:
            return False, f"{len(missing)} entries missing, {len(extra)} extra"
        ref_values = [v for k in ref for v in ref[k]]
        pairs = ((a, b) for k in ref for a, b in zip(new[k], ref[k]))
    limit = atol + rtol * max((abs(v) for v in ref_values), default=0.0)
    diff = max_diff(pairs)
    return diff <= limit, f"max diff {diff:.2e}, limit {limit:.2e}"


def select(data, name, line):
    """output step `line` (counted from 1) of an output, as step 0"""
    if name in DENSE:
        return [data[line - 1]] if len(data) >= line else []
    return {(0, j, i): v for (idx, j, i), v in data.items() if idx == line - 1}


# ----------------------------------------------------------------------
#   running godar
# ----------------------------------------------------------------------

def set_namelist(text, settings):
    """change the value of some namelist entries, keeping the comments"""
    for key, value in settings.items():
        pattern = re.compile(rf"^(\s*{key}\s*=\s*)[^!\n]*?(\s*(!.*)?)$",
                             re.M | re.I)
        text, count = pattern.subn(rf"\g<1>{value}\g<2>", text)
        if count != 1:
            raise ValueError(f"{count} lines set {key} in the namelist")
    return text


def run_godar(args, work, nml, stdin, ranks):
    """run godar in work, return (ok, message, loop time, wall time)"""
    cmd = [str(args.godar)]
    if args.mpiexec:
        cmd = [args.mpiexec, "-n", str(ranks)] + cmd
    elif ranks > 1:
        return False, "--mpiexec is needed for several ranks", None, None
    log = work / f"log.{nml[:-4]}"
    start = time.time()
    with open(log, "w") as out:
        proc = subprocess.run(cmd, input=stdin, text=True, cwd=work,
                              stdout=out, stderr=subprocess.STDOUT,
                              timeout=args.timeout)
    wall = time.time() - start
    text = log.read_text(errors="replace")
    loop = re.search(r"Total simulation time:\s*([0-9.Ee+-]+)", text)
    if proc.returncode != 0 or loop is None:
        return False, f"godar failed, see {log}", None, wall
    # the model falls back on defaults when a namelist cannot be read
    if "error, default will be used" in text:
        return False, f"namelist read error, see {log}", None, wall
    return True, "", float(loop.group(1)), wall


def git_version():
    try:
        rev = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short",
                              "HEAD"], capture_output=True, text=True,
                             check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status",
                                "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True).stdout.strip()
        return rev + ("-dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


# ----------------------------------------------------------------------
#   one test case
# ----------------------------------------------------------------------

def run_case(args, name, case, tolsets):
    """run all the runs of a case, then all its checks. Returns a list of
    (check, ok, messages) and a list of timing rows"""
    work = args.work / name
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(TESTS / "files", work / "files")
    (work / "namelist").mkdir()
    (work / "output").mkdir()

    base = (TESTS / "namelist" / case["namelist"]).read_text()
    expno = {r: f"{k + 1:02d}" for k, r in enumerate(case["runs"])}
    results, timings, failed_runs = [], [], set()

    # runs
    for rname, run in case["runs"].items():
        nml = f"{rname}.nml"
        text = set_namelist(base, run.get("set", {}))
        (work / "namelist" / nml).write_text(text)
        ranks, threads = run.get("ranks", 1), run.get("threads", 1)
        if "restart_from" in run:
            stdin = (f"1\n{nml}\n1\n{expno[run['restart_from']]}\n"
                     f"{run['restart_line']}\n{expno[rname]}\n"
                     f"{case['n']}\n{threads}\n")
        else:
            stdin = f"1\n{nml}\n0\n{expno[rname]}\n{case['n']}\n{threads}\n"

        print(f"  run {rname:8s} ({ranks} rank(s), {threads} thread(s))",
              end="", flush=True)
        ok, msg, loop, wall = run_godar(args, work, nml, stdin, ranks)
        print(f": {loop:.2f} s" if ok else f": FAILED, {msg}")
        if not ok:
            failed_runs.add(rname)
            results.append((f"run {rname}", False, [msg]))
            continue
        steps = re.search(r"^\s*nt\s*=\s*([^!\s]+)", text, re.M | re.I)
        timings.append([name, rname, ranks, threads, case["n"],
                        int(float(steps.group(1).lower().replace("d", "e"))),
                        f"{loop:.3f}", f"{wall:.3f}"])

    def out(rname, fname):
        return work / "output" / f"{fname}.{expno[rname]}"

    # checks
    for check in case["checks"]:
        kind = check["type"]
        runs = check.get("runs", [check.get("run")])
        label = f"{kind} {'/'.join(runs)}"
        if kind == "mass":
            label += f" (max relative change {check['max_rel_change']:g})"
        else:
            label += f" ({check['tolerance']})"
        if failed_runs & set(runs):
            results.append((label, False, ["a run of this check failed"]))
            continue

        msgs, ok_all = [], True
        if kind == "mass":
            r = read_output(out(runs[0], "r"), "r")
            h = read_output(out(runs[0], "h"), "h")
            mass = [sum(a * a * b for a, b in zip(rr, hh))
                    for rr, hh in zip(r, h)]
            change = max(abs(m - mass[0]) for m in mass) / mass[0]
            ok_all = change <= check["max_rel_change"]
            msgs.append(f"relative change {change:.2e}")
            results.append((label, ok_all, msgs))
            continue

        tolset = tolsets[check["tolerance"]]["files"]
        for fname in check.get("files", DENSE + SPARSE):
            tol = tolset.get(fname, tolset.get("*"))
            if tol is None:
                continue
            if kind == "reference":
                ref = TESTS / "reference" / name / f"{fname}.gz"
                if args.update_reference:
                    store_reference(out(runs[0], fname), ref)
                if not ref.exists():
                    ok, msg = False, "no reference, run --update-reference"
                else:
                    ok, msg = compare(read_output(out(runs[0], fname), fname),
                                      read_output(ref, fname), fname, tol)
            elif kind == "same":
                ok, msg = compare(read_output(out(runs[1], fname), fname),
                                  read_output(out(runs[0], fname), fname),
                                  fname, tol)
            elif kind == "restart":
                line = case["runs"][runs[0]]["restart_line"]
                src = case["runs"][runs[0]]["restart_from"]
                ok, msg = compare(
                    select(read_output(out(runs[0], fname), fname), fname, 1),
                    select(read_output(out(src, fname), fname), fname, line),
                    fname, tol)
            else:
                raise ValueError(f"unknown check type {kind}")
            ok_all &= ok
            if args.verbose or not ok:
                msgs.append(f"{fname:7s} {'ok  ' if ok else 'FAIL'} {msg}")
        results.append((label, ok_all, msgs))

    return results, timings


# ----------------------------------------------------------------------
#   main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--godar", type=Path, default=ROOT / "bin" / "godar",
                        help="godar executable (default: bin/godar)")
    parser.add_argument("--mpiexec", default=shutil.which("mpirun")
                        or shutil.which("mpiexec"),
                        help="MPI launcher (default: mpirun from the PATH)")
    parser.add_argument("--work", type=Path, default=TESTS / "work",
                        help="folder for the runs (default: tests/work)")
    parser.add_argument("--case", action="append",
                        help="case to run, can be repeated (default: all)")
    parser.add_argument("--update-reference", action="store_true",
                        help="store the outputs as the new references")
    parser.add_argument("--timeout", type=float, default=3600,
                        help="maximum time of one run [s]")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="show the comparison of every file")
    args = parser.parse_args()
    args.godar = args.godar.resolve()
    args.work = args.work.resolve()

    config = json.loads((TESTS / "cases.json").read_text())
    names = args.case or list(config["cases"])
    version = git_version()
    print(f"godar {args.godar}, version {version}")

    n_failed, n_checks, all_timings = 0, 0, []
    for name in names:
        print(f"\ncase {name}: {config['cases'][name]['description']}")
        results, timings = run_case(args, name, config["cases"][name],
                                    config["tolerances"])
        all_timings += timings
        for label, ok, msgs in results:
            n_checks += 1
            n_failed += not ok
            print(f"  {'PASS' if ok else 'FAIL'}  {label}")
            for msg in msgs:
                print(f"          {msg}")

    # run times, kept from one test session to the next
    log = args.work / "timings.csv"
    new_file = not log.exists()
    with open(log, "a", newline="") as f:
        writer = csv.writer(f)
        if new_file:
            writer.writerow(["date", "version", "host", "case", "run",
                             "ranks", "threads", "n", "steps",
                             "loop_time_s", "wall_time_s"])
        now = datetime.datetime.now().isoformat(timespec="seconds")
        for row in all_timings:
            writer.writerow([now, version, socket.gethostname()] + row)

    if args.update_reference:
        print("\nreferences updated in tests/reference: commit them with "
              "the reason of the change")
    print(f"\n{n_checks - n_failed}/{n_checks} checks passed, "
          f"run times appended to {log}")
    return 1 if n_failed else 0


if __name__ == "__main__":
    sys.exit(main())
