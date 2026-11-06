# Running GODAR on LUMI with amdflang.

After cloning this repository run the scripts under `./lumi`:

* `./build.sh` - build the code depdences as well as the code itself.
* `./run_benchmark.sh` - run an example benchmark that is present under `./tests`.

We currently provide `./build-f08-bug.sh` to replicate an issue arising from mpi_f08 module limitations.
