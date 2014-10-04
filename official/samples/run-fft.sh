#!/bin/bash

EX() {
    echo $@
    $@
}

EX ./fft 1m

echo
echo
EX ulimit -v 50000
EX ./fft 1m

echo
echo
EX ulimit -v 30000
EX ./fft 1m
