#!/bin/sh

./rcdaq2hits /data/gem/xrayscan/COMPASS_4GEM_trial3/COMAPSS_4GEM__0000000108-0000.evt
./hits2centroid COMAPSS_4GEM__0000000108-0000.evt_HITS.root
