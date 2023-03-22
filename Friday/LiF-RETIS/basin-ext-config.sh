#!/bin/bash 
# this script is to extract a configuration right before the lambda 0 interface so that we can use it as a shooting point for the zero minus path 

interface_num="basin"
op_num=1.9 # op for defining the basin boundary 
exttime=`awk -v op_num="$op_num" '{ if ( $2 > op_num ) {print $1} }' init-path.txt | tail -n -1` 
beginLNR=`awk -v op_num="$op_num" '{ if ( $2 > op_num ) {print NR+1} }' init-path.txt | tail -n -1`
begintime=`awk -v beginLNR="$beginLNR" '{ if ( NR == beginLNR ) {print $1; exit} }' init-path.txt`
bop_num=`awk -v beginLNR="$beginLNR" '{ if ( NR == beginLNR ) {print $2; exit} }' init-path.txt`

#echo $op_num $exttime $beginLNR $begintime $bop_num > test.txt

#gmx trjconv -f ~/lj-8192-proj/initial_config/basin/prod.xtc -s ~/lj-8192-proj/initial_config/basin/prod.tpr -o ~/lj-8192-proj/config_extraction/gl-initf_ens-${interface_num}.trr -dump "$exttime" <<< "0"

# because our trajectory is reversed, the -b -e inputs would need to be flipped
gmx trjconv -f init-path.trr -o config_extraction/gl-initf_ens-${interface_num}.trr -b "$exttime" -e "$begintime"

awk -v interface_num="$interface_num" -v begintime="$begintime" -v bop_num="$bop_num" 'BEGIN{print "gl-initf_ens-" interface_num " " begintime " " bop_num " " 1}' >> config_extraction/path_gl-initf_ens-${interface_num}.txt

awk -v interface_num="$interface_num" -v exttime="$exttime" -v op_num="$op_num" 'BEGIN{print "gl-initf_ens-" interface_num " " exttime " " op_num " " 1}' >> config_extraction/path_gl-initf_ens-${interface_num}.txt

