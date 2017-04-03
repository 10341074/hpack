#!/bin/bash

################### TO CHANGE
NAME_FIG_F=fig_ninv
NAME_SRC_F="__mainpb__ __checks__"
#NAME_SRC_2=__checks__
NAME_SRC=runs/src/__mainpb__
NAME_FIG=runs/fig/$NAME_FIG_F
NAME_SRC_DIR=runs/src/

FROM_SRC_DIR=src/
FROM_FIG_DIR=./
#################### NO CHANGE

maxit=0

max(){ # generates sporious files
  if [ $1 > $2 ]
  then
    maxit=$1
  else
    maxit=$2
  fi
}

_OLDS=$(ls $NAME_SRC*)
for i in $_OLDS
do
  it=${i#$NAME_SRC}
  it=${it%.py}
  #max $it $maxit
  if [[ $it -gt $maxit ]] # double [], otherwise creates files?
  then
    maxit=$it
  fi
done

##################### NO CHANGE
maxit=$(($maxit + 1))
echo $maxit
for n in $NAME_SRC_F
do
  #echo $FROM_SRC_DIR''$n''.py $NAME_SRC_DIR''$n''$maxit.py
  cp $FROM_SRC_DIR''$n''.py $NAME_SRC_DIR''$n''$maxit.py
done
## ##################echo $FROM_SRC_DIR''__checks__.py $NAME_SRC''$maxit.py

###################### 

_FIGS=$(ls $FROM_FIG_DIR''$NAME_FIG_F*)

for it in $_FIGS
do
  _EXTS=${it#$FROM_FIG_DIR''$NAME_FIG_F}
   #echo $it $NAME_FIG''$maxit''$_EXTS
   cp $it $NAME_FIG''$maxit''$_EXTS
done


