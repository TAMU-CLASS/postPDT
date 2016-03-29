#! /usr/bin/env bash
# Compare TRIGA pincell PDT to MCNP
#
# tally      description
#  14 -- inner fuel  fis
#  24 -- inner fuel  abs
#  34 -- outer fuel  fis
#  44 -- outer fuel  abs
#  54 -- fuel        fis
#  64 -- fuel        abs

mcnpDir=../dat/mcnp
pdtDir=../dat/pdt

# Run problem
mcnpName=triga_fuel_pin_2d
for prob in {triga_244_136,}
do
    pdtName=out_${prob}
    echo '>>>>' $prob '>>>>'
    ../src/comparePincell.py -m $mcnpName -p $pdtName -i $pdtDir -I $mcnpDir -t 14 34 54 64 -T
done

exit
