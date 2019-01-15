#!/usr/bin/bash

### test to see if the preprocessing code has worked or not
current_dir=$(pwd)
for event_dir in 19* 20*
do
cd "${current_dir}/${event_dir}/raw/"
count=$(ls -1 *.SAC 2>/dev/null | wc -l)

if [ $count != 0 ]; then
  echo "$event_dir didn't work"
else
  :
fi

cd ../../
done
