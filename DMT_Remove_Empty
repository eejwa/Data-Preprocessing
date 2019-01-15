#!/usr/bin/bash

current_dir=$(pwd)

for dir1 in 19* 20* ;do
	cd ${current_dir}/$dir1
	if [ "$(ls -A raw)" ];then
		echo "$dir1 has data, happy days"
		cd ../
	else
		echo "BOOO $dir1 has no data, removing"
		cd ../
		rm -r -f $dir1
	fi
done

## show what's in the file

#for dir2 in 20* ;do
#	cd $dir2
#	if [ "$(ls -A raw)" ];then
##		echo "$dir2 has data, happy days"
#		cd ../
#	else
#		echo "BOOO $dir2 has no data, removing"
#		cd ../
#		rm -r -f $dir2
#	fi
#done
