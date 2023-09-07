#!/bin/bash
iter=0
n_files=$(ls | grep 'F_0*' | wc -l);
nm1_file=$(( $n_files - 1 ))

for f in */ ; do
 iter=$(( $iter + 1 ))
 if [ "$iter" -eq "$nm1_file" ]
 then
  file_get=$f
 fi
done

echo "$file_get" > filename.txt
