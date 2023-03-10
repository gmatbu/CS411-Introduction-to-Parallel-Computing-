#!/bin/bash

damps=('0.10' '0.25' '0.50' '0.90')

# hold k, d constant change d
# d = {.10, .25, .50, .90}
for file in ./data/*.txt; 
do for i in "${damps[@]}";
    do ./main "$file" 8 1000 $i
    done;
done;