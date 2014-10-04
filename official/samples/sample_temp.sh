#!/bin/sh

if [ "$#" -ne 2 ]; then 
    echo "Usage: ./sample.sh mm(=97) fmin(=10 Hz)"
    exit 1
fi

mmin=$1
mmin2=0.$mmin
fmin=$2
fmax=2048 fs=4096 npoint=131072

fn1='gridpoint_hex_all_3pn_vrseb_mm'$mmin'_f'$fmin'.dat'
fn2='vertex_hex_all_3pn_vrseb_mm'$mmin'_f'$fmin'.dat'
fn3='ellipse_hex_all_3pn_vrseb_mm'$mmin'_f'$fmin'.dat'

iDet=2 # vrseb
mm=$mmin2 
m1min=1 m1max=30 m2min=1 m2max=30 m1ini=1.4 m2ini=1.4
nGridMax=250000 nGenMax=10000

./KGLGenerateHexagonalTemplateBank $fn1 $fn2 $fn3 $iDet $fmin $fmax $fs $npoint $mm $m1min $m1max $m2min $m2max $m1ini $m2ini $nGridMax $nGenMax
