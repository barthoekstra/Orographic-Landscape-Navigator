# Orographic-Landscape-Navigator

@NOTE: Add description of the project / abstract

## Requirements
1. MATLAB with:
  1. [Mapping Toolbox](https://www.mathworks.com/products/mapping.html)
  2. [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html). You could do without, but your code will be (much) slower.
2. GDAL command line tools
3. Python3

## Troubleshooting

### 1. GDAL commands not working because of library issues (related to e.g. libtiff)
In order for gdalwarp to work, libtiff from MATLAB may have to be replaced, which can be done by moving the existing library to another folder (libs.bak) and symlinking a new one, as follows:

    cd /Applications/MATLAB_R2016b.app/bin/maci64
    mkdir libs.bak
    mv libtiff.5.dylib libs.bak/
    ln -s /usr/local/opt/libtiff/lib/libtiff.5.dylib

[Source](https://github.com/kyamagu/mexopencv/issues/250#issuecomment-238020778)
