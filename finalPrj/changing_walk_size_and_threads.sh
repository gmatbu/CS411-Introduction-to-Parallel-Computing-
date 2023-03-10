#!/bin/bash

# hold d constant, change k and p
# k = {100, 1000, 10000}
# p = {1, 2, 4, 8}
for i in `seq 2 4`;
    do for j in `seq 0 3`;
    do ./main ./data/facebook_combined.txt $((1 << j)) $((10**$i)) 0.10
    done;
done;