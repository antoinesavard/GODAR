#!/bin/bash

install_name_tool -change @rpath/libcoretran.dylib /opt/coretran/lib/libcoretran.dylib godar

mpirun --bind-to none -n 1 ./godar <input >out
