#!/bin/bash

for argtask in 1 2 3 4 
do

for arg3 in 0.1 0.3 0.5 0.7 0.9
do

for arg2 in 1 2 3 4 6 7
do
if [ "$arg2" = "1" ]; then
	echo Roe, CFL=$arg3
fi
if [ "$arg2" = "2" ]; then
	echo Rusanov, CFL=$arg3
fi
if [ "$arg2" = "3" ]; then
	echo LLF, CFL=$arg3
fi
if [ "$arg2" = "4" ]; then
	echo EO, CFL=$arg3
fi
if [ "$arg2" = "5" ]; then
	echo TVD, CFL=$arg3
fi
if [ "$arg2" = "6" ]; then
	echo BGP, CFL=$arg3
fi
if [ "$arg2" = "7" ]; then
	echo BGP_ch, CFL=$arg3
fi


for arg1 in 32 64 128 256
do
	./main_Roe $argtask $arg1 $arg2 $arg3
done | awk 'BEGIN{enc=0;enl1=0;enl2=0;} {e2nc=$2;e2nl1=$3;e2nl2=$4; pc=log(enc/e2nc)/log(2-1/($1+1)); pl1=log(enl1/e2nl1)/log(2-1/($1+1)); pl2=log(enl2/e2nl2)/log(2-1/($1+1));if(enc!=0){printf("%s %3.3f %3.3f %3.3f\n",$0,pc,pl1,pl2);}else{print($0" - - -");}enc=e2nc; enl1=e2nl1;enl2=e2nl2;}' 


done> ../dat/$argtask/cfl$arg3.dat
done
done