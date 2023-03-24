# this script is used to set up directories that we copy .txt and .trr files into and then used in the global-manager.py

for i in {0..13}; do 
    mkdir ${i}/
    mkdir ${i}/gro
    mkdir ${i}/input
    mkdir ${i}/tpr
    mkdir ${i}/log
    mkdir ${i}/xtc
    mkdir ${i}/trr
    mkdir ${i}/rst
    mkdir ${i}/op
    mkdir ${i}/paths
    cp masters/path_gl-initf_ens-${i}.txt ${i}/paths/path_gl-initf_ens-${i}.txt
    cp masters/gl-initf_ens-${i}*.trr ${i}/trr/
done

