#! /bin/bash


# No need to do this every time
make bernd


# Bernd's script decides on these...
lat_deg1=50.478238579
lon_deg1=8.2900669
lat_deg2=50.47829943178574
lon_deg2=8.289963


# Get original distance
output=`./bernd  --lon_deg1 $lon_deg1 --lat_deg1 $lat_deg1 \
                 --lon_deg2 $lon_deg2 --lat_deg2 $lat_deg2 `


# Sanity check before continuing
#-------------------------------
# Expected outcome is something like
#
# Orig distance [km]: 0.013368735719120646
#
that_worked=1
if [ `echo $output | awk '{print $1}'` != "Orig" ]; then 
    that_worked=0
fi
if [ `echo $output | awk '{print $2}'` != "distance" ]; then 
    that_worked=0
fi
if [ `echo $output | awk '{print $3}'` != "[km]:" ]; then 
    that_worked=0
fi
if [ $that_worked -eq 1 ]; then
    echo "Computation of orig distance worked..."
else
    echo "Computation of orig distance didn't work..."
    exit
fi

command=`echo $output | awk '{print "orig_dist_in_km=\""$4"\""}'`
eval $command

echo "Orig distance: "$orig_dist_in_km



# Now Bernd does his magic and decrees the desired distance:
required_distance_metres=15.0 


# work out coordinates of intermediate point
output=`./bernd  --lon_deg1 $lon_deg1 --lat_deg1 $lat_deg1 \
                  --lon_deg2 $lon_deg2 --lat_deg2 $lat_deg2\
         --required_distance_metres $required_distance_metres` 


# Sanity check before continuing
#-------------------------------
# Expected outcome is something like
#
# Lon (degrees) = 50.47829571941876 Lat (degrees) = 8.2900302707667102 
#

that_worked=1
if [ `echo $output | awk '{print $1}'` != "Lat" ]; then 
    that_worked=0
fi
if [ `echo $output | awk '{print $2}'` != "(degrees)" ]; then 
    that_worked=0
fi
if [ `echo $output | awk '{print $5}'` != "Lon" ]; then 
    that_worked=0
fi
if [ `echo $output | awk '{print $6}'` != "(degrees)" ]; then 
    that_worked=0
fi
if [ $that_worked -eq 1 ]; then
    echo "Computation of new coords worked..."
else
    echo "Computation of new coords didn't work..."
    exit
fi

command=`echo $output | awk '{print "lat_new=\""$4"\" ; lon_new=\""$8"\""}'`
eval $command


# bernd does his magic with these numbers...
echo "Lat of intermediate point: "$lat_new
echo "Lon of intermediate point: "$lon_new

