# this script will make directories for .txt and .trr files to be copied in to it and used in global-manager.py 

mkdir basin/
mkdir basin/gro
mkdir basin/input
mkdir basin/tpr
mkdir basin/log
mkdir basin/xtc
mkdir basin/trr
mkdir basin/rst
mkdir basin/op
mkdir basin/paths
cp masters/path_gl-initf_ens-basin.txt basin/paths/
cp masters/gl-initf_ens-basin*.trr basin/trr/
