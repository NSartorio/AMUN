#!/bin/sh
#
files=`ls $1/*.F90`
for src in $files; do
  deps=`grep -i '  use ' $src | awk '{ print $2 }' | sed 's/,//g' | sort | uniq`
  fname=`basename $src .F90`
  output="$2/${fname}.o: $src"
  for dep in $deps; do
    if [ $dep != "mpi" -a $dep != "hdf5" ]; then
      output="${output} $2/$dep.o"
    fi
  done
  echo $output
done
