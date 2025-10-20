#!/bin/bash

cd ..

mpirun --bind-to none -n 1 ./bin/godar <inputs/input >out
