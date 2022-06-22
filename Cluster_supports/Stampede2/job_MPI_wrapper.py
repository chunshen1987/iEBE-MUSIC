#!/usr/bin/env python3

from mpi4py import MPI
from subprocess import call
import sys
from os import chdir

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

event_id = rank
chdir("event_{}".format(event_id))
call("bash submit_job.script", shell=True)
