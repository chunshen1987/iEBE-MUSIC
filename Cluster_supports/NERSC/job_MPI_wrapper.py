#!/usr/bin/env python3

from mpi4py import MPI
from subprocess import call
import sys
from os import chdir

def print_Usage():
    print("Usage: {} n_threads job_id0 ".format(sys.argv[0]))

try:
    n_threads = int(sys.argv[1])
    job_id    = int(sys.argv[2])
except IndexError:
    print_Usage()
    exit(0)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

event_id = rank + n_threads*(job_id - 1)
chdir("event_{}".format(event_id))
call("bash submit_job.pbs {}".format(rank), shell=True)
