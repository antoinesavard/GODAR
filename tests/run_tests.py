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

Run the tests with ctest (cd build && ctest --output-on-failure), which
calls this script once per case. Only the Python standard library is
needed. To store new references after an intended change of the physics,
from the repository root:

    python tests/run_tests.py --update-reference
"""

import argparse
import csv
import datetime
import gzip
import json
import math
import os
import re
import shutil
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path
from types import SimpleNamespace

TESTS = Path(__file__).resolve().parent
ROOT = TESTS.parent

# outputs with one line of n values per output step
DENSE = ["x", "y", "u", "v", "theta", "omega", "r", "h", "tfx", "tfy",
         "mom", "tsigxx", "tsigyy", "tsigxy", "tsigyx", "tp"]

# outputs with one line "idx j i [values]" per bond or contact
SPARSE = ["bond", "damage", "beam", "angle"]

# what each output holds, printed with the comparisons
QUANTITIES = {
    "x": "position x [m]",
    "y": "position y [m]",
    "u": "velocity x [m/s]",
    "v": "velocity y [m/s]",
    "theta": "angle [rad]",
    "omega": "angular velocity [rad/s]",
    "r": "radius [m]",
    "h": "thickness [m]",
    "tfx": "total force x [N]",
    "tfy": "total force y [N]",
    "mom": "total moment [N m]",
    "tsigxx": "stress xx",
    "tsigyy": "stress yy",
    "tsigxy": "stress xy",
    "tsigyx": "stress yx",
    "tp": "pressure",
    "bond": "list of the bonds",
    "damage": "damage of the bonds [-]",
    "beam": "geometry of the bonds",
    "angle": "contact direction [deg]",
}

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
    """compare two outputs of the same file, new against ref. tol =
    [atol, rtol]: the largest difference over all the values must stay
    below atol + rtol * max|ref|. Returns a namespace: ok, diff (largest
    |new - ref|), scale (max|ref|), limit, count (number of values), and
    problem: why the values cannot be compared one by one, else None"""
    atol, rtol = tol
    result = SimpleNamespace(ok=False, diff=None, scale=None, limit=None,
                             count=0, problem=None)
    if name in DENSE:
        if len(new) != len(ref) or \
                any(len(a) != len(b) for a, b in zip(new, ref)):
            result.problem = (
                f"{len(new)} output steps of {len(new[0]) if new else 0} "
                f"floes, against {len(ref)} of {len(ref[0]) if ref else 0}")
            return result
        ref_values = [v for line in ref for v in line]
        pairs = zip((v for line in new for v in line), ref_values)
    else:
        missing, extra = ref.keys() - new.keys(), new.keys() - ref.keys()
        if missing or extra:
            result.problem = (f"{len(missing)} entries missing, "
                              f"{len(extra)} new entries")
            return result
        ref_values = [v for k in ref for v in ref[k]]
        pairs = ((a, b) for k in ref for a, b in zip(new[k], ref[k]))
    result.count = len(ref_values)
    result.scale = max((abs(v) for v in ref_values), default=0.0)
    result.limit = atol + rtol * result.scale
    result.diff = max_diff(pairs)
    result.ok = result.diff <= result.limit
    return result


def relative(value, scale):
    """value divided by the largest |value| of a file (unchanged when the
    file only holds zeros)"""
    return value / scale if scale else value


def times(ratio):
    """a ratio written short: 23x, 5360x, 3e6x"""
    if ratio < 1e5:
        return f"{ratio:.0f}x"
    return f"{ratio:.0e}x".replace("e+0", "e").replace("e+", "e")


def limit_text(c):
    """limit of a comparison, relative to the largest |value| of the file"""
    lim = relative(c.limit, c.scale)
    return "exact" if lim == 0 else f"{lim:.0e}"


def summarize(comparisons):
    """columns rel. diff, file, limit and margin of a check, for the file
    that came the closest to its limit"""
    problems = [f for f, c in comparisons if c.problem]
    if problems:
        return "-", " ".join(problems), "", "not comparable"

    def closeness(c):
        if c.limit:
            return c.diff / c.limit
        return 0.0 if c.diff == 0 else math.inf

    fname, c = max(comparisons, key=lambda fc: closeness(fc[1]))
    if c.diff == 0:
        exact = all(x.limit == 0 for _, x in comparisons)
        return ("identical", f"{len(comparisons)} files",
                "exact" if exact else "", "")
    rel = f"{relative(c.diff, c.scale):.1e}"
    if c.limit == 0:
        differ = [f for f, x in comparisons if x.diff > 0]
        more = f" +{len(differ) - 1}" if len(differ) > 1 else ""
        return rel, differ[0] + more, "exact", "not identical"
    if c.ok:
        return rel, fname, limit_text(c), times(1 / closeness(c))
    return rel, fname, limit_text(c), times(closeness(c)) + " over"


# one line per run, one line per check, one line per file in the details
RUN_ROW = "  {:9}{:58}{:>9}"
CHECK_ROW = "  {:8}{:25}{:10}{:>9}  {:8}{:>6}  {:>8}"
DETAIL_ROW = "    {:8}{:25}{:>10}{:>10}{:>10}{:>7}  {}"


def detail_line(fname, c):
    """line of a file in the details of a check: largest |value| of the
    file compared against, largest difference, that difference relative
    to the largest |value|, and the limit"""
    quantity = QUANTITIES.get(fname, "")
    if c.problem:
        return f"    {fname:8}{quantity:25}FAIL: {c.problem}"
    verdict = "ok" if c.ok else "FAIL"
    if c.count == 0:
        return (f"    {fname:8}{quantity:25}{'same entries':>30}"
                f"{'':7}  {verdict}")
    return DETAIL_ROW.format(fname, quantity, f"{c.scale:.2e}",
                             f"{c.diff:.2e}",
                             f"{relative(c.diff, c.scale):.1e}",
                             limit_text(c), verdict)


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


def run_godar(args, work, nml, stdin, ranks, cmd=None, env=None):
    """run godar in work, return (ok, message, loop time [s], wall time [s],
    peak memory of the largest process [GB]). cmd: full launch command
    (default: mpiexec -n ranks godar); env: extra environment variables"""
    if cmd is None:
        cmd = [str(args.godar)]
        if args.mpiexec:
            cmd = [args.mpiexec, "-n", str(ranks)] + cmd
        elif ranks > 1:
            return (False, "--mpiexec is needed for several ranks",
                    None, None, None)
    log = work / f"log.{nml[:-4]}"
    start = time.time()
    with open(log, "w") as out:
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True,
                                cwd=work, stdout=out,
                                stderr=subprocess.STDOUT,
                                env=None if env is None
                                else {**os.environ, **env})
        timer = threading.Timer(args.timeout, proc.kill)
        timer.start()
        proc.stdin.write(stdin)
        proc.stdin.close()
        # wait4 also gives the peak resident memory of the largest process,
        # the MPI processes started by mpirun included (bytes on macOS,
        # kB on Linux; on macOS, memory compression makes it vary)
        _, status, usage = os.wait4(proc.pid, 0)
        timer.cancel()
        proc.returncode = os.waitstatus_to_exitcode(status)
    wall = time.time() - start
    memory = usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024) / 1e9
    text = log.read_text(errors="replace")
    loop = re.search(r"Total simulation time:\s*([0-9.Ee+-]+)", text)
    if proc.returncode != 0 or loop is None:
        return False, f"godar failed, see {log}", None, wall, memory
    # the model falls back on defaults when a namelist cannot be read
    if "error, default will be used" in text:
        return False, f"namelist read error, see {log}", None, wall, memory
    return True, "", float(loop.group(1)), wall, memory


def append_rows(path, header, rows):
    """append rows to a csv file; a file with another header (written by
    an older version of the scripts) is kept aside under a dated name"""
    if path.exists():
        with open(path, newline="") as f:
            old = next(csv.reader(f), None)
        if old != header:
            path.rename(path.with_name(
                f"{path.stem}.{datetime.date.today()}{path.suffix}"))
    new_file = not path.exists()
    with open(path, "a", newline="") as f:
        writer = csv.writer(f)
        if new_file:
            writer.writerow(header)
        writer.writerows(rows)


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

def plural(count, word):
    """count and word, e.g. "2 MPI processes", "1 OpenMP thread\""""
    if count == 1:
        return f"{count} {word}"
    return f"{count} {word}{'es' if word.endswith('s') else 's'}"


def run_text(run):
    """MPI processes, OpenMP threads and namelist changes of a run"""
    parts = [plural(run.get("ranks", 1), "MPI process"),
             plural(run.get("threads", 1), "OpenMP thread")]
    parts += [f"{k} = {v}" for k, v in run.get("set", {}).items()]
    return ", ".join(parts)


def differences(a, b):
    """what run b changes with respect to run a"""
    out = []
    for key, word in (("ranks", "MPI process"), ("threads", "OpenMP thread")):
        if a.get(key, 1) != b.get(key, 1):
            out.append(f"{plural(b.get(key, 1), word)} instead of "
                       f"{a.get(key, 1)}")
    sa, sb = a.get("set", {}), b.get("set", {})
    for k in sorted(sa.keys() | sb.keys()):
        if sa.get(k) != sb.get(k):
            out.append(f"{k} = {sb.get(k, 'value of the namelist')}")
    return ", ".join(out)


def run_change(first, run):
    """what a run changes with respect to the first run of its case"""
    if run is first:
        return run_text(run)
    parts = []
    if "restart_from" in run:
        parts.append(f"restarted from output line {run['restart_line']} "
                     f"of {run['restart_from']}")
    if differences(first, run):
        parts.append(differences(first, run))
    return ", ".join(parts) or "nothing"


def comparison(case, check):
    """short name of what a check compares"""
    kind = check["type"]
    if kind == "reference":
        return f"{check['run']} vs stored reference"
    if kind == "same":
        a, b = check["runs"]
        return f"{b} vs {a}"
    if kind == "restart":
        run = case["runs"][check["run"]]
        return (f"{check['run']} vs {run['restart_from']} "
                f"line {run['restart_line']}")
    if kind == "mass":
        return f"mass of {check['run']}"
    raise ValueError(f"unknown check type {kind}")


def run_case(args, name, case, tolsets):
    """run all the runs of a case, printing one line per run, then all
    its checks. Returns the check results (see check_case) and a list of
    timing rows"""
    work = args.work / name
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(TESTS / "files", work / "files")
    (work / "namelist").mkdir()
    (work / "output").mkdir()

    base =(TESTS / "namelist" / case["namelist"]).read_text()
    expno = {r: f"{k + 1:02d}" for k, r in enumerate(case["runs"])}
    results, timings, failed_runs = [], [], set()
    # srun starts the processes on the compute nodes: their memory is not
    # seen from here, only the memory of srun itself
    srun = bool(args.mpiexec) and Path(args.mpiexec).name == "srun"

    print(RUN_ROW.format("run", "what changes", "time")
          + ("" if srun else "    memory"))
    first = next(iter(case["runs"].values()))
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

        print(f"  {rname:9}{run_change(first, run):58}", end="", flush=True)
        ok, msg, loop, wall, memory = run_godar(args, work, nml, stdin,
                                                ranks)
        if not ok:
            print(f"   FAILED, {msg}")
            failed_runs.add(rname)
            results.append(dict(name=f"run {rname}", tolerance=None,
                                ok=False, details=[], note=msg,
                                row=("FAIL", f"run {rname}", "", "", "",
                                     "", "failed")))
            continue
        print(f"{loop:7.2f} s" + ("" if srun else f"   {memory:.2f} GB"))
        steps = re.search(r"^\s*nt\s*=\s*([^!\s]+)", text, re.M | re.I)
        timings.append([name, rname, ranks, threads, case["n"],
                        int(float(steps.group(1).lower().replace("d", "e"))),
                        f"{loop:.3f}", f"{wall:.3f}",
                        "" if srun else f"{memory:.3f}"])

    results += check_case(args, name, case, tolsets, work, expno,
                          failed_runs)
    return results, timings


def check_case(args, name, case, tolsets, work, expno, failed_runs):
    """all the checks of a case, on the outputs of its runs in work.
    Returns one dict per check: name (what is compared), tolerance, ok,
    row (columns of the table of the checks), details (the comparison of
    each file) and note (why the check could not be made)"""
    def out(rname, fname):
        return work / "output" / f"{fname}.{expno[rname]}"

    results = []
    for check in case["checks"]:
        kind = check["type"]
        runs = check.get("runs", [check.get("run")])
        what = comparison(case, check)
        tname = check.get("tolerance")
        result = dict(name=what, tolerance=tname, ok=False, details=[],
                      note=None, row=None)
        results.append(result)

        if failed_runs & set(runs):
            result["note"] = "not checked: a run of this check failed"
            result["row"] = ("FAIL", what, tname or "-", "", "", "",
                             "not run")
            continue

        if kind == "mass":
            r = read_output(out(runs[0], "r"), "r")
            h = read_output(out(runs[0], "h"), "h")
            mass = [sum(a * a * b for a, b in zip(rr, hh))
                    for rr, hh in zip(r, h)]
            change = max(abs(m - mass[0]) for m in mass) / mass[0]
            limit = check["max_rel_change"]
            result["ok"] = change <= limit
            margin = "" if change == 0 else (
                times(limit / change) if result["ok"]
                else times(change / limit) + " over")
            result["row"] = ("PASS" if result["ok"] else "FAIL", what, "-",
                             f"{change:.1e}", "r**2 h", f"{limit:.0e}",
                             margin)
            continue

        tolset = tolsets[tname]["files"]
        comparisons = []
        for fname in check.get("files", DENSE + SPARSE):
            tol = tolset.get(fname, tolset.get("*"))
            if tol is None:
                continue
            if kind == "reference":
                ref = TESTS / "reference" / name / f"{fname}.gz"
                if args.update_reference:
                    store_reference(out(runs[0], fname), ref)
                if not ref.exists():
                    c = SimpleNamespace(
                        ok=False, problem="no reference, run "
                        "--update-reference")
                else:
                    c = compare(read_output(out(runs[0], fname), fname),
                                read_output(ref, fname), fname, tol)
            elif kind == "same":
                c = compare(read_output(out(runs[1], fname), fname),
                            read_output(out(runs[0], fname), fname),
                            fname, tol)
            elif kind == "restart":
                line = case["runs"][runs[0]]["restart_line"]
                src = case["runs"][runs[0]]["restart_from"]
                c = compare(
                    select(read_output(out(runs[0], fname), fname), fname, 1),
                    select(read_output(out(src, fname), fname), fname, line),
                    fname, tol)
            else:
                raise ValueError(f"unknown check type {kind}")
            comparisons.append((fname, c))

        result["ok"] = all(c.ok for _, c in comparisons)
        result["details"] = comparisons
        result["row"] = ("PASS" if result["ok"] else "FAIL", what, tname,
                         *summarize(comparisons))
    return results


def print_checks(results, verbose):
    """table of the checks of a case, then the details of the failed
    checks (of all of them when verbose)"""
    print(CHECK_ROW.format("result", "comparison", "tolerance", "rel. diff",
                           "file", "limit", "margin"))
    for result in results:
        print(CHECK_ROW.format(*result["row"]).rstrip())
    for result in results:
        if result["note"]:
            print(f"\n  {result['name']}: {result['note']}")
        if result["details"] and (verbose or not result["ok"]):
            print(f"\n  details of {result['name']} "
                  f"(tolerance {result['tolerance']}):")
            print(DETAIL_ROW.format("file", "quantity", "max|value|",
                                    "max diff", "rel. diff", "limit",
                                    "").rstrip())
            for fname, c in result["details"]:
                print(detail_line(fname, c))


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
    parser.add_argument("-v", "-V", "--verbose", action="store_true",
                        help="list every file in the details of every "
                        "check, not only of the failed ones")
    args = parser.parse_args()
    args.godar = args.godar.resolve()
    args.work = args.work.resolve()

    config = json.loads((TESTS / "cases.json").read_text())
    names = args.case or list(config["cases"])
    version = git_version()
    print(f"godar {args.godar}, version {version}")

    failed, n_checks, all_timings = [], 0, []
    for name in names:
        case = config["cases"][name]
        print(f"\ncase {name}: {case['description']}\n")
        results, timings = run_case(args, name, case, config["tolerances"])
        all_timings += timings
        print()
        print_checks(results, args.verbose)
        for result in results:
            n_checks += 1
            if not result["ok"]:
                failed.append(f"{name} {result['name']}")

    # run times, kept from one test session to the next
    log = args.work / "timings.csv"
    now = datetime.datetime.now().isoformat(timespec="seconds")
    append_rows(log, ["date", "version", "host", "case", "run", "ranks",
                      "threads", "n", "steps", "loop_time_s", "wall_time_s",
                      "peak_memory_GB"],
                [[now, version, socket.gethostname()] + row
                 for row in all_timings])

    if args.update_reference:
        print("\nreferences updated in tests/reference: commit them with "
              "the reason of the change")
    print(f"\n{n_checks - len(failed)}/{n_checks} checks passed, "
          f"run times appended to {log}")
    if failed:
        print("failed: " + ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
