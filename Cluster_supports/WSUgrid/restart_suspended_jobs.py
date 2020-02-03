#!/usr/bin/env python3

import sys
from os import path
import subprocess
from datetime import datetime, timedelta

try:
    jobid = str(sys.argv[1])
except:
    print("Usage: {} job_id".format(str(sys.argv[0])))
    exit(1)

cmd = subprocess.Popen(['qstat', '-f', '-w', jobid], stdout=subprocess.PIPE)
cmd_out, cmd_err = cmd.communicate()
job_info = cmd_out.decode(sys.stdout.encoding).split("\n")
#print(job_info)

job_state = ''
start_date = ''
job_path = ''
for line_i in job_info:
    if "job_state" in line_i:
        job_state = line_i.split("=")[1].strip()
        if job_state != "S":
            print("Job is not suspended, job state = {}".format(job_state))
            print("exit")
            exit(0)
    if "etime" in line_i:
        start_date = line_i.split("=")[1].strip()
        start_date = datetime.strptime(start_date, '%a %b %d %H:%M:%S %Y')
    if "Output_Path" in line_i:
        job_path_list = line_i.split("=")[1].strip().split(":")[1].split("/")
        job_path = "/".join(job_path_list[:-1])

check_date = datetime.today()
print("check_date : ",check_date)
print("etime : ", start_date)
print("job path: ", job_path)

if (check_date - start_date > timedelta(days=1)):
    print("restart job: ", jobid)
    cmd = subprocess.Popen(['qdel', jobid], stdout=subprocess.PIPE)
    subprocess.call("(cd {}; qsub submit_job.pbs;)".format(job_path), shell=True)
