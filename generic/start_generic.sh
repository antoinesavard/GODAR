#!/bin/bash

mpirun --bind-to none -n 1 ./godar <input_restart >out
