% S = shaperead('data/ahn/ahn_units.shp');
% mapshow(S)

profile on
info = geotiffinfo('data/ahn/r14bn2.tif')
[X, R] = geotiffread('data/ahn/r14bn2.tif');
%%
Z = resizem(X, 0.2);

%%
R_Half = R;
R_half.RasterSize = size(Z);

geotiffwrite('data/ahn/r14bn2b.tif', Z, R_Half, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
%imwrite(Z,'data/ahn/r14bn2b.tif');
profile viewer

% Command used:
% gdalwarp r14bn2.tif r14bn2c.tif -t_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" -ts 2000 2500  -r bilinear

% Changed to:
% gdalwarp r14bn2a.tif r14bn2a2.tif -t_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" -ts 2000 2500  -r bilinear

% Start over:
% First test the resizing visually:
% gdalwarp r14bn2.tif r14bn2a-small.tif -ts 2000 2500 -r bilinear

% Then reproject
% gdalwarp r14bn2a-small.tif r14bn2a-small-rp.tif -t_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

%%
% Now all at once
% Clearly the best option
% Really
% gdalwarp r14bn2.tif r14bn2-small-rp-2.tif -ts 2000 2500 -r bilinear -t_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

%% Let's try reprojecting the shapefile
% Doesn't work:
% ogr2ogr -f "ESRI Shapefile" ahn_units_wgs84.shp ahn_units.shp -t_srs EPSG:4326

% Important to place output file before input file
% ogr2ogr -s_srs EPSG:28992 -t_srs EPSG:4326 ahn_units_wgs84.shp ahn_units.shp

%% Clip ahn_units_wgs84 by area selected
% ogr2ogr -clipsrc ResearchAreaNH2.shp out.shp ahn_units_wgs84.shp
% setenv('PATH', [getenv('PATH') ':/usr/local/bin']); % Change to real PATH
setenv('PATH', '/usr/local/bin:/Users/barthoekstra/anaconda/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/Library/TeX/texbin');

Selection = shaperead('data/out.shp');
Units = string({Selection.UNIT});
% http://geodata.nationaalgeoregister.nl/ahn2/extract/ahn2_05m_ruw/r09ez1.tif.zip
base_url = 'http://geodata.nationaalgeoregister.nl/ahn2/extract/ahn2_05m_ruw/r%s.tif.zip'
base_filename = 'data/dem/ahn2_05m_ruw_%s.tif.zip'
for i = 1:2 %numel(Units)
    url = sprintf(base_url,Units(i))
    filename = sprintf(base_filename,Units(i))
    outfilename = websave(filename,url) % Need to create dem folder first? Yes
    outfilename2 = unzip(outfilename,'data/dem/')
    delete(char(outfilename));
    base_reproject = 'gdalwarp "%s" "data/dem/r%s.wgs84.tif" -ts 2000 2500 -r bilinear -t_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"'
    cmd = sprintf(base_reproject, char(outfilename2), Units(i))
    cd(pwd)
    [status,cmdout] = system(cmd)
    delete(char(outfilename2))
end

% In order for gdalwarp to work, libtiff from MATLAB may have to be
% replaced, which can be done by moving the existing library and symlinking
% a new one, as follows:
%
% cd /Applications/MATLAB_R2016a.app/bin/maci64
% mkdir libs.bak
% mv libtiff.5.dylib libs.bak/
% ln -s /usr/local/opt/libtiff/lib/libtiff.5.dylib
% 
% Source: https://github.com/kyamagu/mexopencv/issues/250#issuecomment-238020778

%% Shapereads
shapeinfo('data/ResearchAreaNH2.shp')
S = shaperead('data/ResearchAreaNH2.shp')

shapeinfo('data/ahn/ahn_units.shp')
S2 = shaperead('data/ahn/ahn_units.shp')

shapeinfo('data/ahn/ahn_units_wgs84.shp')
S3 = shaperead('data/ahn/ahn_units_wgs84.shp')

%% Selection
SXY = [S.X;S.Y]'

arr = [];
for i = 1:size(S3,1)
    for j = 1:numel(SXY)
        if SXY(j) == 
    end
end
