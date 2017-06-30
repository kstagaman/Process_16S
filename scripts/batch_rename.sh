#!/bin/bash
# batch_rename.sh

if [ $# != 3 ]
then
	echo "Usage: batch_rename.sh <file identifier> <part to rename> <substitution>"
else
	for file in $1
	do 
		mv $file `echo $file | sed "s/$2/$3/"`
	done 
fi
