# Custom environment shell code should follow

# Disabled - please use osgvo-el7-cuda10 if you want this functionality
#if [ "x$LD_LIBRARY_PATH" = "x" ]; then
#    export LD_LIBRARY_PATH="/host-libs"
#else
#    export LD_LIBRARY_PATH="/host-libs:$LD_LIBRARY_PATH"
#fi

# ensure we have PS1 set
PS1="Singularity $SINGULARITY_NAME:\\w> "
export PS1

