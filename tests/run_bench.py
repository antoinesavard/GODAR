#!/usr/bin/env python3
"""
GODAR benchmarks.

Times godar on the cases of tests/cases.json, with longer runs than the
tests and outputs only at the start and the end. With --baseline, two
executables are compared (A/B mode): their runs alternate (A B, B A,
A B, ...) so that both see the same state of the machine, and the ratio
of their times is computed for each pair of runs, which is much less
sensitive to the machine than comparing times measured at different
moments. Without --baseline, only --godar is timed.

The time of each phase of a time step comes from the timers printed by
godar (pkg/timers); executables without them only give the total time.
The results are appended to <work>/bench.csv.

Only the Python standard library is needed. From the repository root:

    python tests/run_bench.py --baseline /path/to/other/godar
    python tests/run_bench.py --case TESTSNAP --steps 500 --repeat 5
    python tests/run_bench.py --threads 4
"""

import argparse
import csv
import datetime
import json
import re
import shutil
import socket
import statistics
from pathlib import Path
from types import SimpleNamespace

import run_tests as rt


def time_run(exe, args, work, nml, stdin):
    """run one executable, return (loop time, {phase: time})"""
    run_args = SimpleNamespace(godar=exe, mpiexec=args.mpiexec,
                               timeout=args.timeout)
    ok, msg, loop, _ = rt.run_godar(run_args, work, nml, stdin, args.ranks)
    if not ok:
        raise SystemExit(f"{exe}: {msg}")
    log = (work / f"log.{nml[:-4]}").read_text(errors="replace")
    phases = {m.group(1): float(m.group(2)) for m in
              re.finditer(r"^\s*Timer\s+(\w+)\s*:\s*([0-9.]+)", log, re.M)}
    return loop, phases


def spread(values):
    """median (min-max) of a list of values"""
    return (f"{statistics.median(values):8.3f} "
            f"({min(values):.3f}-{max(values):.3f})")


def report(results, labels):
    """table of the median time of each phase, and of the ratios
    candidate / baseline of the pairs of runs"""
    phases = ["total"] + list(results["candidate"][0][1])
    header = f"  {'phase':8s}" + "".join(f"{l:>26s}" for l in labels)
    if "baseline" in labels:
        header += f"{'candidate/baseline':>26s}"
    print(header)
    for phase in phases:
        values = {}
        for label in labels:
            if phase == "total":
                values[label] = [loop for loop, _ in results[label]]
            elif all(phase in p for _, p in results[label]):
                values[label] = [p[phase] for _, p in results[label]]
        line = f"  {phase:8s}"
        for label in labels:
            line += f"{spread(values[label]) if label in values else '-':>26s}"
        if "baseline" in values and "candidate" in values:
            ratios = [c / b for c, b in zip(values["candidate"],
                                             values["baseline"]) if b > 0]
            if ratios:
                line += f"{spread(ratios):>26s}"
        print(line)

    if "baseline" in labels:
        ratios = [c / b for (c, _), (b, _) in zip(results["candidate"],
                                                  results["baseline"])]
        med = statistics.median(ratios)
        if min(ratios) <= 1 <= max(ratios):
            print("  -> no difference beyond the run-to-run variation "
                  f"(paired ratios {min(ratios):.3f}-{max(ratios):.3f})")
        else:
            word = "faster" if med < 1 else "slower"
            print(f"  -> candidate {abs(1 - med) * 100:.1f}% {word} "
                  f"(paired ratios {min(ratios):.3f}-{max(ratios):.3f})")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--godar", type=Path,
                        default=rt.ROOT / "bin" / "godar",
                        help="executable to time (default: bin/godar)")
    parser.add_argument("--baseline", type=Path,
                        help="executable to compare with (A/B mode)")
    parser.add_argument("--mpiexec", default=shutil.which("mpirun")
                        or shutil.which("mpiexec"),
                        help="MPI launcher (default: mpirun from the PATH)")
    parser.add_argument("--case", action="append",
                        help="case to run, can be repeated (default: all)")
    parser.add_argument("--steps", type=int, default=1000,
                        help="time steps of each run (default: 1000)")
    parser.add_argument("--repeat", type=int, default=3,
                        help="runs of each executable (default: 3)")
    parser.add_argument("--ranks", type=int, default=1,
                        help="MPI processes (default: 1)")
    parser.add_argument("--threads", type=int, default=1,
                        help="OpenMP threads (default: 1)")
    parser.add_argument("--work", type=Path, default=rt.TESTS / "work"
                        / "bench", help="folder for the runs "
                        "(default: tests/work/bench)")
    parser.add_argument("--timeout", type=float, default=7200,
                        help="maximum time of one run [s]")
    args = parser.parse_args()

    exes = {"candidate": args.godar.resolve()}
    if args.baseline:
        exes["baseline"] = args.baseline.resolve()
    labels = list(exes)
    expno = {label: f"{k + 1:02d}" for k, label in enumerate(labels)}

    config = json.loads((rt.TESTS / "cases.json").read_text())
    version = rt.git_version()
    print(f"candidate {exes['candidate']} (repository version {version})")
    if args.baseline:
        print(f"baseline  {exes['baseline']}")

    log = args.work / "bench.csv"
    for name in args.case or list(config["cases"]):
        case = config["cases"][name]
        print(f"\ncase {name}: {args.steps} steps, {args.ranks} rank(s) x "
              f"{args.threads} thread(s), {args.repeat} runs of each")

        work = args.work / name
        if work.exists():
            shutil.rmtree(work)
        shutil.copytree(rt.TESTS / "files", work / "files")
        (work / "namelist").mkdir()
        (work / "output").mkdir()
        text = rt.set_namelist((rt.TESTS / "namelist" / case["namelist"])
                               .read_text(), {"nt": f"{args.steps}d0",
                                              "comp": f"{args.steps}d0"})
        for label in labels:
            (work / "namelist" / f"{label}.nml").write_text(text)

        results = {label: [] for label in labels}
        for rep in range(args.repeat):
            for label in labels if rep % 2 == 0 else labels[::-1]:
                stdin = (f"1\n{label}.nml\n0\n{expno[label]}\n{case['n']}\n"
                         f"{args.threads}\n")
                loop, phases = time_run(exes[label], args, work,
                                        f"{label}.nml", stdin)
                results[label].append((loop, phases))
                print(f"  run {rep + 1} {label:9s}: {loop:.3f} s")
        report(results, labels)

        new_file = not log.exists()
        with open(log, "a", newline="") as f:
            writer = csv.writer(f)
            if new_file:
                writer.writerow(["date", "version", "host", "case", "steps",
                                 "ranks", "threads", "label", "executable",
                                 "run", "loop_time_s", "phases_s"])
            now = datetime.datetime.now().isoformat(timespec="seconds")
            for label in labels:
                for rep, (loop, phases) in enumerate(results[label]):
                    writer.writerow([now, version, socket.gethostname(),
                                     name, args.steps, args.ranks,
                                     args.threads, label, exes[label],
                                     rep + 1, f"{loop:.3f}",
                                     ";".join(f"{k}={v:.3f}"
                                              for k, v in phases.items())])
    print(f"\nresults appended to {log}")


if __name__ == "__main__":
    main()
