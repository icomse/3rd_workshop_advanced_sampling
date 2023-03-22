#!/bin/bash

MYPATH=${1}
MPATH=${2}

mkdir ${MPATH}
mkdir ${MPATH}data
mkdir ${MPATH}simulations

cp ${MYPATH}data/*.itp ${MPATH}data
cp ${MYPATH}data/*.ndx ${MPATH}data
cp ${MYPATH}data/*.top ${MPATH}data
cp ${MYPATH}data/*.mdp ${MPATH}data
