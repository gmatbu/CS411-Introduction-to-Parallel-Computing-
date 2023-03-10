#!/bin/bash

# hold w, d constant change p = {1 , 2, 4, 8}
for file in ./data/*.txt; 
do for i in `seq 0 3`;
    do ./main "$file" $((1 << i)) 1000 0.10
    done;
done;