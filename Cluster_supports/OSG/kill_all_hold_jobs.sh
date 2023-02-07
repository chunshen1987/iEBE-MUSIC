#!/usr/bin/env bash

condor_q -hold | grep "chunshen1987 " | awk {'print $1'} | xargs condor_rm
