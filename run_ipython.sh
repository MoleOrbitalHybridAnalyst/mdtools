#!/bin/bash

#module="python/2.7-2015q2"
module="python"

if [ "$#" -eq 1 ]; then
   module=$1
fi

module load $module
module load spark/2.1
export PYTHONPATH=$SPARK_HOME/python:$PYTHONPATH

#ip=$(host `uname -n` | cut -d ' ' -f 4)
ip=$(/sbin/ip route get 8.8.8.8 | awk '{print $NF;exit}')
port=$((10000+ $RANDOM % 20000))

echo "Loaded Python and Spark modules"
echo "Starting ipython notebook ..."


jupyter notebook --no-browser --ip=$ip --port=$port --log-level='ERROR'
