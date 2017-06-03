#!/usr/bin/python

# This code is simply a wrapper for running gdal commands, without MATLAB
# causing issues with dependencies, etc.

import sys
import os

print(sys.argv[0])

action = sys.argv[1]
targetfile = sys.argv[2]

if action == "merge":
    print('Mergeing...')
    # gdalbuildvrt merged.vrt r14bn2.wgs84.tif r14en1.wgs84.tif r14ez1.wgs84.tif r14bz2.wgs84.tif r14bz1.wgs84.tif r14bn1.wgs84.tif r09dz1.wgs84.tif r09dz2.wgs84.tif r09gz1.wgs84.tif
    # gdalbuildvrt output.vrt files-to-be-merged.tif separated-by-spaces.tif
    # python gisactions.py merge data/dem/output.vrt data/dem/r14bn2.wgs84.tif data/dem/r14bn1.wgs84.tif

    # First create a virtual mosaic
    # Cmd format: gdalbuildvrt output.vrt file1.tif file2.tif file3.tif
    print('Creating mosaic...')
    targetvrt = targetfile.replace(".tif", ".vrt")
    cmd_mosaic = "gdalbuildvrt %s %s" % (targetvrt, ' '.join(sys.argv[3:]))
    os.system(cmd_mosaic)

    # Now translate the mosaic to an actual GeoTiff
    # Cmd format: gdal_translate -of GTiff mosaic.vrt output.tif
    mergedfile = sys.argv[2].replace(".wgs84.tif", ".merged.wgs84.vrt")
    cmd_merge = "gdal_translate -a_srs \"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs\" -of GTiff %s %s" % (targetvrt, targetfile)
    os.system(cmd_merge)

    # Now remove the .vrt
    os.remove(targetvrt)

    print('Merge finished...')
elif action == "reproject":
    print('Reprojecting...')
else:
    print('No valid action provided.')
