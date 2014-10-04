#!/bin/sh

fn1=test_f256-512.dat
fn2=sample.ana
fs=1024 tmin=0 tmax=10
nframe=256 nstart=0 nend=512 nshift=200 Fthr=-8 a2thr=0.01 nsig=20 nitr=1000 mu0=1 nu0=1

./nha $fn1 $fn2 $fs $tmin $tmax $nframe $nstart $nend $nshift $offset $Fthr $a2thr $nsig $nitr $mu0 $nu0
