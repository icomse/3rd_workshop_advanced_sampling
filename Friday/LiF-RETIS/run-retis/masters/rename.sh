# this script is to rename the files so that it's compatible with global-manager.py's function called read_paths

for i in {0..9}
do
    mv gl-initf_ens-0${i}.trr gl-initf_ens-${i}.trr
    mv path_gl-initf_ens-0${i}.txt path_gl-initf_ens-${i}.txt
    oldname="gl-initf_ens-0${i}"
    newname="gl-initf_ens-${i}"
    echo $oldname
    echo $newname
    cp path_gl-initf_ens-${i}.txt path_gl-initf_ens-${i}.txt.bak
    sed "s/$oldname/$newname/g" path_gl-initf_ens-${i}.txt.bak > path_gl-initf_ens-${i}.txt
done


