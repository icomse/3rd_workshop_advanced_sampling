#!/bin/bash 

mkdir config_extraction 

# now actually finding the time of the op value

while read line; 
do
    interface_num=`awk '{print $1}' <<< "$line"`
    op_num=`awk '{print $2}' <<< "$line"`
    exttime=`awk -v op_num="$op_num" '{ if ( $2 > op_num ) {print $1} }' init-path.txt | tail -n -1`

    gmx trjconv -f init-path.trr -o config_extraction/gl-initf_ens-${interface_num}.trr -dump "$exttime" <<< "0" 
    
    awk -v interface_num="$interface_num" -v exttime="$exttime" -v op_num="$op_num" 'BEGIN{print "gl-initf_ens-" interface_num " " exttime " " op_num " " 1}' >> config_extraction/path_gl-initf_ens-${interface_num}.txt <<< "$line"

    echo "gl-initf_ens-${interface_num} 0.0 0 1" >> config_extraction/path_gl-initf_ens-${interface_num}.txt 

done < trial_interfaces.txt

