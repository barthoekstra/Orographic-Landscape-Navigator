# Orographic-Landscape-Navigator

## Thesis abstract
Orographic lift, the upward deflection of horizontal airflows, can be caused by large-scale structures, such as mountain ridges, but also by micro-topographic features, such as dunes, dykes, buildings and even vegetation. Gulls, due to their morphology, are able to continuously adjust their flight strategy to the local aerial conditions. The extent to which these flight generalists use small-scale orographic lift conditions to optimise their energy efficiency of flight, has not been studied before. Research has shown, however, that a significant share of soaring occurs at low altitudes, suggesting a strong link to surface processes. The aim of this study was to quantify the opportunities and conditions for soaring on orographic lift along the flyways of a Dutch population of herring gulls (Larus argentatus) and lesser black-backed gulls (L. fuscus). Flight tracks were annotated with orographic lift rates, modelled using a high-resolution digital surface model. Low-altitude soaring is preferred when orographic lift rates are between 0.8 and 4.2 m/s. In herring gulls a larger share of this soaring can be explained by orographic lift than in lesser black-backed gulls, 35% and 19% respectively, which indicates a stronger preference for orographic lift soaring in the former species. Predicting low-altitude soaring using logistic regression models suggests a modest influence of orographic lift on the soaring flight probability, as the best model was only up to 60% accurate. This, however, is most likely a result of the quality of the used predictors, rather than an indication of the limited use of orographic lift to facilitate low- altitude soaring.

## Requirements
1. MATLAB with:
  1. [Mapping Toolbox](https://www.mathworks.com/products/mapping.html)
  2. [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html). You could do without, but your code will be (much) slower, in which case you need to replace all occurences of `parfor` with `for`.
2. GDAL command line tools
3. Python3

## Data sources
1. Digital Surface Model from the [Algemeen Hoogtebestand Nederland](http://www.ahn.nl/index.html)
2. Hourly weather (wind) data from the [Dutch Meteorological Institute (KNMI)](http://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script)
3. Weather station data from [data.overheid.nl](https://data.overheid.nl/data/dataset/a4c85871-6847-4b55-bc46-e14d74cf2c05)

All data acquisition is automated and thoroughly documented within the code.

## Troubleshooting

### 1. GDAL commands not working because of library issues (related to e.g. libtiff)
In order for gdalwarp to work, libtiff from MATLAB may have to be replaced, which can be done by moving the existing library to another folder (libs.bak) and symlinking a new one, as follows:

    cd /Applications/MATLAB_R2016b.app/bin/maci64
    mkdir libs.bak
    mv libtiff.5.dylib libs.bak/
    ln -s /usr/local/opt/libtiff/lib/libtiff.5.dylib

[Source](https://github.com/kyamagu/mexopencv/issues/250#issuecomment-238020778)

If the solution above does still not solve issues, you could rewrite the code to generate a simple text file with all the commands listed. It is then possible to just copy-paste this list of commands in the Terminal and execute the statements. This is inconvenient, but so is fixing issues related to the interaction of MATLAB and GDAL libraries.
