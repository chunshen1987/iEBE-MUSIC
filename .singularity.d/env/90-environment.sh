# Custom environment shell code should follow

if [ "x$LD_LIBRARY_PATH" = "x" ]; then
    export LD_LIBRARY_PATH="/host-libs"
else
    export LD_LIBRARY_PATH="/host-libs:$LD_LIBRARY_PATH"
fi

