#!/bin/bash

source env.config

VCLANG=${VAN}/clang
MCLANG=${MOD}/clang

HOME=`pwd`/../NPB3.3-SER-C

VBIN=`pwd`/vanilla_bins
MBIN=`pwd`/auto2_bins

BLD=`pwd`/build

mkdir -p $VBIN $MBIN $BLD

CFLAGS="-Oz -fno-crash-diagnostics"
LFLAGS=

TARGETS="bt cg ep ft is lu mg sp"
CLASS=S

# configure extension file here for testing
EXTENSION=

printf "========= BUILDING VANILLA TARGETS ==========================================\n"
for t in $TARGETS; do
    printf "\n\nBuilding auto2 binaries for $t and class $CLASS\n"
    ./build.sh $t $CLASS $HOME $MBIN $BLD $VCLANG $MCLANG "$CFLAGS" "$LFLAGS" "$EXTENSION"
done


printf "\n\n========= BUILDING AUTO2 TARGETS ============================================\n"
EXTENSION=""

for t in $TARGETS; do
    printf "\n\nBuilding vanilla binaries for $t and class $CLASS\n"
    ./build.sh $t $CLASS $HOME $VBIN $BLD $VCLANG $VCLANG "$CFLAGS" "$LFLAGS" "$EXTENSION"
done