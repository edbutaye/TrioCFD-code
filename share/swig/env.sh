#!/bin/bash

# MEDCoupling charge par defaut est l'optim, possible de charger $TRUST_MEDCOUPLING_ROOT/env_debug.sh
. $TRUST_MEDCOUPLING_ROOT/env.sh

export PYTHONPATH=$project_directory/share/swig/install/lib:$PYTHONPATH
export LD_LIBRARY_PATH=$project_directory/build/src/exec_opt:$LD_LIBRARY_PATH  # Not needed anymore ?

