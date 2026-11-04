#!/usr/bin/env bash
# Collects the results of iEBE-MUSIC runs on Roihu (CSC). The collection runs
# in the container, because the login nodes have no Python with h5py.
#
# Run it in the folder where the jobs were generated (the campaign folder for
# a campaign), which holds collect_events_singularity.sh and
# combine_multiple_hdf5.py:
#
#     ./collect_results.sh <image> <run folder> [<run folder> ...]
#
# The results of the event folders of each run folder (RUN_NAME of a job
# array, or a design point of a campaign) are added to
# RESULTS/<run folder>/<run folder>.h5, and the parameter files of the run
# folder are copied next to it. The results of an event folder are moved
# into this file, so collecting again adds only the event folders that have
# finished since. If the collection fails, the results are kept in
# RESULTS/<run folder>/SPVN_RESULTS, and collecting again retries.

usage="usage: $0 <image> <run folder> [<run folder> ...]"

image=$1
shift
if [ -z "$image" ] || [ $# -eq 0 ]; then
    echo "$usage"
    exit 1
fi
if [ ! -f "$image" ]; then
    echo "image $image not found"
    exit 1
fi
if [ ! -f collect_events_singularity.sh ]; then
    echo "run $0 in the folder where the jobs were generated"
    exit 1
fi

status=0
for runFolder in "$@"; do
    runFolder=${runFolder%/}
    runName=`basename $runFolder`
    nFolders=`ls -d $runFolder/event_* 2> /dev/null | wc -l`
    nResults=`ls $runFolder/event_*/RESULTS_*.h5 2> /dev/null | wc -l`
    if [ "$nFolders" -eq 0 ]; then
        echo "$runFolder: no event folders"
        status=1
        continue
    fi
    # results of an earlier collection that failed
    kept=RESULTS/$runName/SPVN_RESULTS
    if [ "$nResults" -eq 0 ] && [ ! -d $kept ]; then
        echo "$runFolder: no new results"
        continue
    fi
    if [ -d $kept ]; then
        echo "$runFolder: collecting $nResults of $nFolders event folders" \
             "and the results kept from a failed collection"
    else
        echo "$runFolder: collecting $nResults of $nFolders event folders"
    fi
    if ! apptainer exec --bind="/scratch,/projappl" "$image" \
            ./collect_events_singularity.sh $runFolder RESULTS; then
        echo "$runFolder: collection failed, no results were deleted;" \
             "collect again to retry"
        status=1
        continue
    fi
    # the parameter files of the run, next to its results
    find $runFolder -maxdepth 1 -type f ! -name "*.script" \
        -exec cp {} RESULTS/$runName/ \;
    apptainer exec --bind="/scratch,/projappl" "$image" python3 -c \
        'import sys, h5py; print(sys.argv[1] + ":", len(h5py.File(sys.argv[1], "r")), "events")' \
        RESULTS/$runName/$runName.h5 || status=1
done
exit $status
