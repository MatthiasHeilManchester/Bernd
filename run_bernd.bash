#! /bin/bash


# No need to do this every time
make bernd


# Bernd's script decides no these...
lon_deg1=50.478238579
lat_deg1=8.2900669
lon_deg2=50.47829943178574
lat_deg2=8.289963


# Get original distance
command=`./bernd  --lon_deg1 $lon_deg1 --lat_deg1 $lat_deg1 \
                  --lon_deg2 $lon_deg2 --lat_deg2 $lat_deg2 | awk '{print "orig_dist_in_km=\""$4"\""}'` 
eval $command

echo "Orig distance: "$orig_dist_in_km

# Now Bernd does his magic and decrees the desired distance:
required_distance_metres=12.0 


# work out coordinates of intermediate point
command=`./bernd  --lon_deg1 $lon_deg1 --lat_deg1 $lat_deg1 \
                  --lon_deg2 $lon_deg2 --lat_deg2 $lat_deg2\
         --required_distance_metres $required_distance_metres | awk '{print "lon_new=\""$4"\" ; lat_new=\""$8"\""}'` 
eval $command


echo "Lon of intermediate point: "$lon_new
echo "Lat of intermediate point: "$lat_new

