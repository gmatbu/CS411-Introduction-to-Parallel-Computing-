#!/bin/bash

# hold p, d constant change k = {100, 1000, 10000}
for file in ./data/*.txt; 
do for i in `seq 2 4`;
    do ./main "$file" 8 $((10**$i)) 0.10
    done;
done;