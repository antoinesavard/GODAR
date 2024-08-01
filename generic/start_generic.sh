#!/bin/bash

cd ..

mpirun --bind-to none -n 1 ./godar <inputs/input >out
