# Running GODAR on LUMI with amdflang

After cloning this repository run the scripts under `./lumi`:

* `./build.sh` - build the code depdences as well as the code itself.

Then, from the folder above `src/`, submit the jobs that run ctest:

* `sbatch src/lumi/ctest_tests.sh` - the regression tests, on the debug partition.
* `sbatch src/lumi/ctest_bench.sh` - the ARCTIC benchmark (`ctest -C bench -R arctic`), on 1 to 4 CPU nodes of the standard partition.

We currently provide `./build-f08-bug.sh` to replicate an issue arising from mpi_f08 module limitations.
