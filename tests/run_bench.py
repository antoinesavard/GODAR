#!/usr/bin/env python3
"""
GODAR benchmarks.

Times godar on the cases of tests/cases.json with longer runs than the
tests (outputs only at the start and the end), for configurations of
nodes x MPI processes per node x OpenMP threads per process. For each
configuration it gives the median time, the speedup and the parallel
efficiency compared with the first configuration, the peak memory of the
largest process, the load imbalance between the processes (slowest /
fastest force loop), and the time of each phase of a time step (timers
printed by godar, pkg/timers).

Larger problems: --particles N1,N2,... runs the ARCTIC case, N floes
filling the NSIDC sea ice extent of the Arctic, touching each other. The
case of each N is built once by tests/bench/make_arctic.py and kept in
<work>/cases. The sizes of the floes, the
search radius and the time step scale with N, so that each floe has the
same neighbours whatever N: the work of a step grows like N (weak
scaling: N times the floes on N times the cores should take the same
time, compare the lines with the same particles per core).

Clusters: inside a Slurm allocation (salloc or sbatch), each
configuration is launched with srun and the threads are pinned to cores
(OMP_PLACES=cores, OMP_PROC_BIND=close). --launch replaces the launch
command for clusters that need something else. The memory needed by
each configuration (about 100 n^2 bytes per process) is checked against
the memory of a node, and configurations that do not fit are skipped.
--dry-run prints the commands without running them, --slurm-script
writes a job script running the same benchmark.

With --baseline, two executables are compared (A/B mode): their runs
alternate (A B, B A, ...) so that both see the same state of the
machine, and the ratio of their times is computed for each pair of runs.
--save-baseline stores the median times as the baseline of the machine
(tests/bench/baselines/<machine>.json); later sessions on the same
machine are compared with it and warn when slower than --threshold %.

The results are appended to <work>/bench.csv. Only the Python standard
library is needed. From the repository root:

    python tests/run_bench.py --threads 1,2,4
    python tests/run_bench.py --baseline /path/to/other/godar
    salloc --nodes=8 --exclusive ...    then
    python tests/run_bench.py --particles 20000 \\
        --nodes 1,2,4,8 --ranks-per-node 1,2 --threads 32,64
"""

import argparse
import datetime
import json
import os
import re
import shlex
import shutil
import socket
import statistics
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import run_tests as rt

BASELINES = rt.TESTS / "bench" / "baselines"
MAKE_ARCTIC = rt.TESTS / "bench" / "make_arctic.py"

# launch commands, the {fields} are filled for each configuration and the
# executable is added at the end
SRUN = (
    "srun --nodes={nodes} --ntasks={ranks} "
    "--ntasks-per-node={ranks_per_node} --cpus-per-task={threads} "
    "--cpu-bind=cores"
)
MPIRUN = "{mpiexec} -n {ranks}"

# options only used to write the Slurm script, not passed to the job
SLURM_OPTIONS = [
    "--slurm-script",
    "--account",
    "--partition",
    "--time",
    "--setup",
]

IN_SLURM = "SLURM_JOB_ID" in os.environ


# ----------------------------------------------------------------------
#   machine and configurations
# ----------------------------------------------------------------------


def int_list(text):
    return [int(v) for v in text.split(",")]


def memory_estimate(n, threads):
    """memory of one process [GB]: the n x n pair state (about 100 n**2
    bytes), the per-particle arrays and their copies for each thread"""
    return (100 * n**2 + (2000 + 120 * threads) * n) / 1e9


def node_memory(args):
    """memory of a node [GB]: --node-memory, else the memory Slurm gives
    the job on a node, else the memory of this machine"""
    if args.node_memory:
        return args.node_memory
    slurm = os.environ.get("SLURM_MEM_PER_NODE")  # MiB
    if slurm and int(slurm) > 0:
        return int(slurm) * 1.048576e-3
    if IN_SLURM:
        return None  # e.g. --mem=0: unknown here, may be a login node
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9
    except (ValueError, OSError):
        return None


def available(args):
    """(nodes, cores per node) available to the benchmark"""
    if IN_SLURM:
        nodes = int(os.environ.get("SLURM_JOB_NUM_NODES", "1"))
        cpus = os.environ.get("SLURM_CPUS_ON_NODE")
        return nodes, int(cpus) if cpus else None
    return 1, os.cpu_count()


def skip_reason(args, config, n, memory, nodes_avail, cpus_avail):
    """why a configuration cannot run, None if it can"""
    nodes, rpn, threads = config
    if nodes > nodes_avail:
        if not IN_SLURM and not args.launch:
            return "several nodes need a Slurm allocation or --launch"
        return f"{nodes} nodes, {nodes_avail} available"
    if cpus_avail and rpn * threads > cpus_avail:
        return f"{rpn * threads} cores per node, {cpus_avail} available"
    need = rpn * memory_estimate(n, threads)
    if memory and need > 0.9 * memory:
        return f"needs about {need:.0f} GB per node, {memory:.0f} GB"
    return None


def launch_command(args, exe, config):
    """command running exe in a configuration (nodes, ranks/node, threads)"""
    nodes, rpn, threads = config
    template = args.launch or (SRUN if IN_SLURM else MPIRUN)
    fields = dict(
        nodes=nodes,
        ranks=nodes * rpn,
        ranks_per_node=rpn,
        threads=threads,
        mpiexec=args.mpiexec or "mpirun",
    )
    return shlex.split(template.format(**fields)) + [str(exe)]


def thread_env(args, threads):
    """OpenMP environment of a configuration"""
    env = {"OMP_NUM_THREADS": str(threads)}
    if args.bind != "none":
        env.update(OMP_PLACES="cores", OMP_PROC_BIND=args.bind)
    return env


# ----------------------------------------------------------------------
#   cases
# ----------------------------------------------------------------------


def arctic_dir(args, requested):
    return args.work / "cases" / f"ARCTIC_{requested}"


def arctic_case(args, requested, build=True):
    """ARCTIC case with about `requested` floes, built if missing (unless
    build is False); returns its folder and its case.json, None if not
    built"""
    folder = arctic_dir(args, requested)
    if not (folder / "case.json").exists():
        if not build:
            return folder, None
        print(f"\nbuilding the ARCTIC case with {requested} floes in {folder}")
        subprocess.run(
            [
                sys.executable,
                str(MAKE_ARCTIC),
                "--n",
                str(requested),
                "--out",
                str(folder),
            ],
            check=True,
        )
    return folder, json.loads((folder / "case.json").read_text())


def jobs_of(args, config):
    """(name, size, files folder, namelist, particles, built) of each case
    to run"""
    jobs = []
    if args.case or not args.particles:
        for name in args.case or list(config["cases"]):
            case = config["cases"][name]
            jobs.append(
                (
                    name,
                    "default",
                    rt.TESTS / "files",
                    rt.TESTS / "namelist" / case["namelist"],
                    case["n"],
                    True,
                )
            )
    for requested in args.particles or []:
        folder, info = arctic_case(args, requested, build=not args.dry_run)
        jobs.append(
            (
                "ARCTIC",
                f"n{requested}",
                folder / "files",
                folder / "namelist" / "ARCTIC.nml",
                info["n"] if info else requested,
                info is not None,
            )
        )
    return jobs


# ----------------------------------------------------------------------
#   runs
# ----------------------------------------------------------------------


def time_run(exe, args, work, nml, stdin, config):
    """run one executable, return a dict with the loop time, the time of
    each phase (rank 0), the peak memory and the load imbalance"""
    nodes, rpn, threads = config
    run_args = SimpleNamespace(
        godar=exe, mpiexec=args.mpiexec, timeout=args.timeout
    )
    ok, msg, loop, _, memory = rt.run_godar(
        run_args,
        work,
        nml,
        stdin,
        nodes * rpn,
        cmd=launch_command(args, exe, config),
        env=thread_env(args, threads),
    )
    if not ok:
        raise SystemExit(f"{exe}: {msg}")
    log = (work / f"log.{nml[:-4]}").read_text(errors="replace")
    timers = re.findall(
        r"^\s*Timer\s+(\w+)\s*:\s*([0-9.]+)\s*\(\s*[0-9.]+"
        r"\s*%\)\s*([0-9.]+)\s+([0-9.]+)",
        log,
        re.M,
    )
    phases = {name: float(t) for name, t, _, _ in timers}
    imbalance = None
    for name, _, tmin, tmax in timers:
        if name == "pairs" and float(tmin) > 0:
            imbalance = float(tmax) / float(tmin)
    return dict(loop=loop, phases=phases, memory=memory, imbalance=imbalance)


def label_of(config):
    return "x".join(str(v) for v in config)


def spread(values, width=24):
    """median (min-max) of a list of values"""
    text = (
        f"{statistics.median(values):.3f} "
        f"({min(values):.3f}-{max(values):.3f})"
    )
    return f"{text:>{width}s}"


# ----------------------------------------------------------------------
#   reports
# ----------------------------------------------------------------------


def summary(results, configs, n, stored, threshold):
    """one line per configuration; returns the warnings"""
    warnings = []
    first = configs[0]
    cores_first = first[0] * first[1] * first[2]
    t_first = statistics.median(r["loop"] for r in results["candidate"][first])
    header = (
        f"  {'NxRxT':10s}{'cores':>6s}{'part/core':>10s}"
        f"{'time (s)':>24s}{'speedup':>8s}{'effic.':>7s}"
        f"{'memory':>9s}{'imbal.':>7s}"
    )
    if "baseline" in results:
        header += f"{'baseline (s)':>24s}{'cand./base.':>24s}"
    if stored:
        header += f"{'stored (s)':>11s}{'ratio':>7s}"
    print(header)

    for config in configs:
        runs = results["candidate"][config]
        cores = config[0] * config[1] * config[2]
        times = [r["loop"] for r in runs]
        speedup = t_first / statistics.median(times)
        efficiency = speedup * cores_first / cores
        imbal = [r["imbalance"] for r in runs if r["imbalance"]]
        line = (
            f"  {label_of(config):10s}{cores:6d}{n / cores:10.0f}"
            f"{spread(times)}{speedup:8.2f}{efficiency * 100:6.0f}%"
            f"{max(r['memory'] for r in runs):6.1f} GB"
            + (f"{statistics.median(imbal):7.2f}" if imbal else f"{'-':>7s}")
        )
        if "baseline" in results:
            base = [r["loop"] for r in results["baseline"][config]]
            line += spread(base) + spread([c / b for c, b in zip(times, base)])
        if stored:
            ref = stored.get(label_of(config), {}).get("total")
            if ref:
                ratio = statistics.median(times) / ref
                line += f"{ref:11.3f}{ratio:7.3f}"
                if ratio > 1 + threshold / 100:
                    warnings.append(
                        f"{label_of(config)}: "
                        f"{(ratio - 1) * 100:.0f}% slower than "
                        "the stored baseline"
                    )
        print(line)
    return warnings


def phase_table(results, config):
    """median time of each phase (rank 0), and the paired ratios
    candidate/baseline"""
    labels = list(results)
    phases = ["total"] + list(results["candidate"][config][0]["phases"])
    print(f"\n  phases, {label_of(config)} (s, rank 0)")
    print(
        f"  {'phase':8s}"
        + "".join(f"{l:>24s}" for l in labels)
        + (f"{'cand./base.':>24s}" if "baseline" in labels else "")
    )
    for phase in phases:
        values = {}
        for label in labels:
            runs = results[label][config]
            if phase == "total":
                values[label] = [r["loop"] for r in runs]
            elif all(phase in r["phases"] for r in runs):
                values[label] = [r["phases"][phase] for r in runs]
        line = f"  {phase:8s}"
        for label in labels:
            line += spread(values[label]) if label in values else f"{'-':>24s}"
        if "baseline" in values and "candidate" in values:
            ratios = [
                c / b
                for c, b in zip(values["candidate"], values["baseline"])
                if b > 0
            ]
            if ratios:
                line += spread(ratios)
        print(line)


def write_slurm_script(args):
    """job script running this benchmark (without the Slurm options) on
    the largest number of nodes asked"""
    argv, skip = [], False
    for a in sys.argv[1:]:
        if skip:
            skip = False
        elif a in SLURM_OPTIONS:
            skip = True
        elif not any(a.startswith(o + "=") for o in SLURM_OPTIONS):
            argv.append(a)
    lines = [
        "#!/bin/bash",
        "#SBATCH --job-name=godar-bench",
        f"#SBATCH --nodes={max(args.nodes)}",
        "#SBATCH --exclusive",
        "#SBATCH --mem=0",
        f"#SBATCH --time={args.time}",
        f"#SBATCH --output={args.work}/slurm-%j.out",
    ]
    if args.account:
        lines.append(f"#SBATCH --account={args.account}")
    if args.partition:
        lines.append(f"#SBATCH --partition={args.partition}")
    lines += ["", "# environment of the cluster (modules, ...)"]
    lines.append(f"export TOPDIR={shlex.quote(str(rt.ROOT.parent))}")
    if args.setup:
        lines += Path(args.setup).read_text().splitlines()
    lines += [
        "",
        f"cd {shlex.quote(str(rt.ROOT))}",
        shlex.join(["python3", str(Path(__file__).resolve())] + argv),
    ]
    Path(args.slurm_script).write_text("\n".join(lines) + "\n")
    args.work.mkdir(parents=True, exist_ok=True)  # Slurm does not create it
    print(
        f"job script written in {args.slurm_script}: "
        f"sbatch {args.slurm_script}"
    )
    missing = [
        n
        for n in args.particles or []
        if arctic_case(args, n, build=False)[1] is None
    ]
    if missing:
        # the masks need ncdump, scipy, numba and the NSIDC download,
        # which compute nodes may not have
        print(
            "the job would first build these ARCTIC cases, better build "
            "them before:"
        )
        for n in missing:
            print(
                "    "
                + shlex.join(
                    [
                        sys.executable,
                        str(MAKE_ARCTIC),
                        "--n",
                        str(n),
                        "--out",
                        str(arctic_dir(args, n)),
                    ]
                )
            )


# ----------------------------------------------------------------------
#   main
# ----------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--godar",
        type=Path,
        default=rt.ROOT / "bin" / "godar",
        help="executable to time (default: bin/godar)",
    )
    parser.add_argument(
        "--baseline", type=Path, help="executable to compare with (A/B mode)"
    )
    parser.add_argument(
        "--case",
        action="append",
        help="test case to run, can be repeated (default: "
        "all, none with --particles)",
    )
    parser.add_argument(
        "--particles",
        type=int_list,
        help="ARCTIC case with about these numbers of floes, e.g. 5000,20000",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=1000,
        help="time steps of each run (default: 1000)",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=3,
        help="runs of each configuration (default: 3)",
    )
    parser.add_argument(
        "--nodes",
        type=int_list,
        default=[1],
        help="nodes, e.g. 1,2,4,8 (default: 1)",
    )
    parser.add_argument(
        "--ranks-per-node",
        "--ranks",
        type=int_list,
        default=[1],
        help="MPI processes per node, e.g. 1,2 (default: 1)",
    )
    parser.add_argument(
        "--threads",
        type=int_list,
        default=[1],
        help="OpenMP threads per process, e.g. 1,2,4,8 (default: 1)",
    )
    parser.add_argument(
        "--launch",
        help="launch command with the fields "
        "{nodes} {ranks} {ranks_per_node} {threads} "
        "{mpiexec}, the executable is added at the end "
        f"(default: '{SRUN}' in Slurm, '{MPIRUN}' outside)",
    )
    parser.add_argument(
        "--mpiexec",
        default=shutil.which("mpirun") or shutil.which("mpiexec"),
        help="MPI launcher outside Slurm (default: mpirun)",
    )
    parser.add_argument(
        "--bind",
        default="close" if IN_SLURM else "none",
        choices=["close", "spread", "none"],
        help="OpenMP thread binding (default: close in Slurm, none outside)",
    )
    parser.add_argument(
        "--node-memory",
        type=float,
        help="memory of a node [GB], to skip configurations "
        "that do not fit (default: from Slurm or this "
        "machine)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the commands without running them",
    )
    parser.add_argument(
        "--save-baseline",
        action="store_true",
        help="store the times as the baseline of the machine",
    )
    parser.add_argument(
        "--machine",
        default=socket.gethostname(),
        help="name of the machine for the stored baseline "
        "(default: host name)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=10,
        help="slowdown compared with the stored baseline "
        "that gives a warning [%%] (default: 10)",
    )
    parser.add_argument(
        "--work",
        type=Path,
        default=rt.TESTS / "work" / "bench",
        help="folder for the runs (default: tests/work/bench)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=7200,
        help="maximum time of one run [s]",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="phase table of every configuration",
    )
    slurm = parser.add_argument_group("job script (--slurm-script)")
    slurm.add_argument(
        "--slurm-script",
        help="write a Slurm job script "
        "running this benchmark, instead of running it",
    )
    slurm.add_argument("--account", help="Slurm account")
    slurm.add_argument("--partition", help="Slurm partition")
    slurm.add_argument(
        "--time",
        default="1:00:00",
        help="time limit of the job (default: 1:00:00)",
    )
    slurm.add_argument(
        "--setup",
        help="file with the lines setting up the "
        "environment in the job (module load ...)",
    )
    args = parser.parse_args()
    args.work = args.work.resolve()

    if args.slurm_script:
        write_slurm_script(args)
        return

    exes = {"candidate": args.godar.resolve()}
    if args.baseline:
        exes["baseline"] = args.baseline.resolve()
    labels = list(exes)
    expno = {label: f"{k + 1:02d}" for k, label in enumerate(labels)}
    configs = [
        (nd, r, t)
        for nd in args.nodes
        for r in args.ranks_per_node
        for t in args.threads
    ]

    config = json.loads((rt.TESTS / "cases.json").read_text())
    version = rt.git_version()
    nodes_avail, cpus_avail = available(args)
    memory = node_memory(args)
    print(f"candidate {exes['candidate']} (repository version {version})")
    if args.baseline:
        print(f"baseline  {exes['baseline']}")
    print(
        f"{'Slurm job ' + os.environ['SLURM_JOB_ID'] if IN_SLURM else 'no Slurm'}"
        f": {nodes_avail} node(s), {cpus_avail or '?'} cores and "
        f"{f'{memory:.0f} GB' if memory else '? GB'} per node, "
        f"thread binding {args.bind}"
    )

    stored_file = BASELINES / f"{args.machine}.json"
    stored = (
        json.loads(stored_file.read_text())
        if stored_file.exists()
        else {"cases": {}}
    )
    if stored_file.exists() and not args.save_baseline:
        print(
            f"stored baseline of {args.machine}: version "
            f"{stored.get('version')}, {stored.get('date')}"
        )
    if args.save_baseline and stored.get("steps") != args.steps:
        # times of another number of steps cannot be compared: start over
        stored = {"cases": {}}

    rows, warnings = [], []
    now = datetime.datetime.now().isoformat(timespec="seconds")
    for name, size, files, namelist, n, built in jobs_of(args, config):
        work = args.work / (name if size == "default" else f"{name}_{size}")
        print(
            f"\ncase {name}: {n} particles"
            f"{'' if built else ' (the case is not built yet)'}, "
            f"{args.steps} steps, {args.repeat} runs of each "
            "configuration (nodes x ranks/node x threads)"
        )

        todo = []
        for c in configs:
            why = skip_reason(args, c, n, memory, nodes_avail, cpus_avail)
            if why:
                print(f"  skip {label_of(c)}: {why}")
            else:
                todo.append(c)
        if args.dry_run:
            for c in todo:
                print(
                    f"  {label_of(c)}: about "
                    f"{c[1] * memory_estimate(n, c[2]):.1f} GB per "
                    f"node, {' '.join(launch_command(args, exes['candidate'], c))}"
                    f"  ({' '.join(f'{k}={v}' for k, v in thread_env(args, c[2]).items())})"
                )
            continue
        if not todo:
            continue

        if work.exists():
            shutil.rmtree(work)
        shutil.copytree(files, work / "files", symlinks=True)
        (work / "namelist").mkdir()
        (work / "output").mkdir()
        text = rt.set_namelist(
            namelist.read_text(),
            {"nt": f"{args.steps}d0", "comp": f"{args.steps}d0"},
        )
        for label in labels:
            (work / "namelist" / f"{label}.nml").write_text(text)

        results = {label: {c: [] for c in todo} for label in labels}
        for c in todo:
            for rep in range(args.repeat):
                for label in labels if rep % 2 == 0 else labels[::-1]:
                    stdin = f"1\n{label}.nml\n0\n{expno[label]}\n{n}\n{c[2]}\n"
                    res = time_run(
                        exes[label], args, work, f"{label}.nml", stdin, c
                    )
                    results[label][c].append(res)
                    print(
                        f"  {label_of(c)} run {rep + 1} {label:9s}: "
                        f"{res['loop']:.3f} s, {res['memory']:.2f} GB"
                    )
                    rows.append(
                        [
                            now,
                            version,
                            socket.gethostname(),
                            name,
                            size,
                            n,
                            args.steps,
                            *c,
                            label,
                            exes[label],
                            rep + 1,
                            f"{res['loop']:.3f}",
                            f"{res['memory']:.3f}",
                            f"{res['imbalance']:.3f}"
                            if res["imbalance"]
                            else "",
                            ";".join(
                                f"{k}={v:.3f}"
                                for k, v in res["phases"].items()
                            ),
                        ]
                    )

        print()
        same = stored.get("steps") == args.steps and not args.save_baseline
        ref = stored["cases"].get(name, {}).get(size, {}) if same else {}
        warnings += [
            f"{name} {size} {w}"
            for w in summary(results, todo, n, ref, args.threshold)
        ]
        for c in todo if args.verbose else todo[:1]:
            phase_table(results, c)

        if args.save_baseline:
            entry = stored["cases"].setdefault(name, {}).setdefault(size, {})
            for c, runs in results["candidate"].items():
                entry[label_of(c)] = {
                    "total": statistics.median(r["loop"] for r in runs),
                    "phases": {
                        k: statistics.median(r["phases"][k] for r in runs)
                        for k in runs[0]["phases"]
                    },
                    "memory_GB": max(r["memory"] for r in runs),
                }

    if args.dry_run:
        return
    log = args.work / "bench.csv"
    rt.append_rows(
        log,
        [
            "date",
            "version",
            "host",
            "case",
            "size",
            "n",
            "steps",
            "nodes",
            "ranks_per_node",
            "threads",
            "label",
            "executable",
            "run",
            "loop_time_s",
            "peak_memory_GB",
            "imbalance",
            "phases_s",
        ],
        rows,
    )

    if args.save_baseline:
        stored.update(
            {
                "machine": args.machine,
                "version": version,
                "date": now,
                "steps": args.steps,
            }
        )
        BASELINES.mkdir(parents=True, exist_ok=True)
        stored_file.write_text(json.dumps(stored, indent=4) + "\n")
        print(f"\nbaseline of {args.machine} stored in {stored_file}")
    for w in warnings:
        print(f"WARNING: {w}")
    print(f"\nresults appended to {log}")


if __name__ == "__main__":
    main()
