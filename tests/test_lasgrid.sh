#!/bin/bash 

# Tests for lasgrid.

txt2las -a_srs epsg:2956 -i data/stats_coords.txt -o /tmp/stats_coords.las

echo "Min"
lasgrid -v -o /tmp/lasgrid_min.tif -t min -r 1 -d -1 -s 2956 /tmp/stats_coords.las

echo "Max"
lasgrid -v -o /tmp/lasgrid_max.tif -t max -r 1 -d -1 -s 2956 /tmp/stats_coords.las

echo "Mean"
lasgrid -v -o /tmp/lasgrid_mean.tif -t mean -r 1 -d -1 -s 2956 /tmp/stats_coords.las

echo "Median"
lasgrid -o /tmp/lasgrid_median.tif -t median -r 1 -d -1 -s 2956 /tmp/stats_coords.las

echo "Variance"
lasgrid -v -o /tmp/lasgrid_variance.tif -t variance -r 1 -d -1 -s 2956 /tmp/stats_coords.las

echo "Pop. Variance"
lasgrid -v -o /tmp/lasgrid_pvariance.tif -t pvariance -r 1 -d -1 -s 2956 /tmp/stats_coords.las

echo "Std Dev"
lasgrid -v -o /tmp/lasgrid_stddev.tif -t stddev -r 1 -d -1 -s 2956 /tmp/stats_coords.las

echo "Pop. Std Dev"
lasgrid -v -o /tmp/lasgrid_pstddev.tif -t pstddev -r 1 -d -1 -s 2956 /tmp/stats_coords.las

echo "Density"
lasgrid -v -o /tmp/lasgrid_density.tif -t density -r 1 -d -1 -s 2956 /tmp/stats_coords.las

echo "Count"
lasgrid -v -o /tmp/lasgrid_count.tif -t count -r 1 -d -1 -s 2956 /tmp/stats_coords.las

echo "Done"
