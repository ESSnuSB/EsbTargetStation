#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64

##################
# ROOT
##################
export ROOTSYS=/libcern/root/5.30.04/sl5.5-x86_64
export PATH=/$ROOTSYS/bin:$PATH:/home/baussan/Tools:.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/$ROOTSYS/lib


source env_493_64_libcern.sh
./SPL_SUPERBEAM run.mac

