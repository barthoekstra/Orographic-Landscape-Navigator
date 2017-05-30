function [tracksOroLift] = calcOrographicLift(storpath, filename, trackselection, stations, wind, cellsize, radius, nanthreshold)
% storpath = 'data/dem/'; 
% filename = 'r09az2.wgs84.tif';
% deminfo = dems(3,:);
% trackselection = tracks(tracks.longitude >= deminfo.XMin & tracks.longitude <= deminfo.XMax ...
%                             & tracks.latitude >= deminfo.YMin & tracks.latitude <= deminfo.YMax, :);
% cellsize = 2.5;
% radius = 4; 
% nanthreshold = -1e10;
%calcOrographicLift Calculates orographic lift on provided locations
%   Inputs:
%   1. Path where DEM tiles are stored (e.g. 'data/dem/')
%   2. Cellsize in meters
%   3. Radius around latlon coordinate to select maximum orographic lift
%      from in cells (e.g. with cellsize of 2.5m, a radius of 4 corresponds
%      with 10m)
%   4. 

%
%   setenv('PATH', ['/usr/local/bin:/Users/barthoekstra/anaconda/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/Library/TeX/texbin:', getenv('PATH')])

% Some global variables
%nanthreshold = -1e10;
%cellsize = 2.5; % m (meter)
%radius = 4; % 10m (meter)
%storpath = 'data/dem/';
filesuffix = '.wgs84.tif';
filename = char(filename);

% ---------------------------------------------------------------------
% Merge AHN tiles so orographic lift can also be calculated along edges
% ---------------------------------------------------------------------

% Calculating the orographic lift in the middle of the DEM tile is simple,
% but it gets more complicated along the edges where possibly no data is
% available (due to the projection) or we need to find a value from a
% selection of raster points which falls outside of the DEM tile. The only
% solution for this problem is to merge the tile with all surrounding tiles
% (when available) and use these to fill in the datagaps
tiles = struct2table(shaperead('data/ahn_units_wgs84.clipped.shp'));
ref_tile_name = replace(filename,filesuffix,'');
ref_tile_name = ref_tile_name(2:end);

ref_tile = tiles(strcmp(ref_tile_name, tiles.UNIT),:);

merge_units = [];

for i = 1:size(tiles,1)
   comp_tile = tiles(i,:); % comparison tile
   
   % Intersect polygons
   [xi, yi] = polyxpoly(cell2mat(ref_tile.X), cell2mat(ref_tile.Y), ...
                        cell2mat(comp_tile.X), cell2mat(comp_tile.Y));
   
   if ~isempty(xi) & ~isempty(yi)
       % Apparently there is some intersection. Store the name of this
       % unit, so we know what to merge later on
       merge_units = [merge_units, comp_tile.UNIT];
   end
   
end

% Now merge the units we have found above as follows:
% 1. Make a virtual mosaic of selected (intersecting) files
% 2. Turn virtual mosaic into a GeoTiff
%
% Luckily, calling gisactions.py with the following command does that all
% at once (hopefully without causing issues with MATLAB dependencies/libs)
files = join(strcat(storpath, 'r', merge_units, filesuffix), ' ');
output = strcat(storpath, 'r', ref_tile_name, '.merged', filesuffix);
cmd = sprintf('python gisactions.py merge %s %s', output, files);
          
[status, cmd] = system(cmd);
if status == 1
    error(cmd)
end

original_dem = [storpath, filename];
merged_dem = output;

%%
% ---------------------------------------------------------------------
% Select all datapoints within boundary of the original (non-merged) 
% DEM
% ---------------------------------------------------------------------
% info = deminfo(strcmp(deminfo.filenames, filename),:);
% 
% trackselection = tracks(tracks.longitude >= info.XMin & tracks.longitude <= info.XMax ...
%                         & tracks.latitude >= info.YMin & tracks.latitude <= info.YMax, :);

%%
% ---------------------------------------------------------------------
% Calculate orographic lift
% ---------------------------------------------------------------------
dem = geotiffread(merged_dem);
dem(dem <= nanthreshold) = 0;

[r, c] = size(dem);

% dzdy
dzdy = zeros(r, c);
dzdy(1,:) = (dem(1,:) - dem(2,:)) ./ (cellsize);
dzdy(r,:) = (dem(r-1,:) - dem(r,:)) ./ (cellsize);
dzdy(2:r-1,:) = (dem(1:r-2,:) - dem(3:r,:)) ./ (cellsize * 2);

% dzdx
dzdx = zeros(r, c);
dzdx(:,1) = (dem(:,1) - dem(:,2)) ./ (cellsize);
dzdx(:,c) = (dem(:,c-1) - dem(:,c)) ./ (cellsize);
dzdx(:,2:c-1) = (dem(:,1:c-2) - dem(:,3:c)) ./ (cellsize * 2);

% Slope
slope = atand(sqrt(dzdx.^2 + dzdy.^2));

% Aspect
aspect = 90 + atan2d(dzdy, dzdx);
aspect = aspect + (aspect < 0) * 360; % transform to [0 360]

% Cleanup
slope(isnan(slope)) = 0;

% @NOTE: Remove the commented out section below later if not used
% Calculate Ca (Updraft coefficient) for every 10 degrees change of wind
% direction. We don't need to calculate Ca in finer steps, because KNMI
% only reports the values in steps of 10 degrees.
% Ca = zeros(size(dem, 1), size(dem, 2), 36);
% for i = 10:10:360
%     k = i / 10;
%     Ca(:,:,k) = sind(slope) .* cosd(i - aspect);
% end

%%
% Now it is - finally - time to calculate the orographic lift
spatialinfo = geotiffinfo(merged_dem);
n = size(trackselection, 1);

% Make some variables to store column data in, so we can merge them with
% the table later on
wspeed_col   = zeros(n,1);
wdir_col     = wspeed_col;
oroglift_col = wspeed_col;
wstation_col = wspeed_col;

parfor i = 1:n
    % Change date formatting to get wind data
    datetime = datenum(trackselection(i,:).date_time, 'yyyy-mm-dd HH:MM:SS');
    date = str2num(datestr(datetime, 'yyyymmdd'));
    
    tracksel = trackselection(i,:);
    lat = tracksel.latitude;
    lon = tracksel.longitude;
    hour = tracksel.hour + 1;
    
    cWS = getNearestWeatherStation(stations, lat, lon); % closest weather station
    
    windselection = wind(wind.stationID == cWS & wind.date == date ...
                         & wind.hour == hour, :);
    wspeed = windselection.wspeed_hr / 10; % convert 0.1 m/s to m/s
    wdir   = windselection.wdir;
    
    [r, c] = setpostn(dem, spatialinfo.SpatialRef, lat, lon);
    rows = r-radius:1:r+radius;
    rows(rows <= 0 | rows > size(dem, 1)) = []; % Cannot select out of range of DEM, so remove these rows
    cols = c-radius:1:c+radius;
    cols(cols <= 0 | cols > size(dem, 2)) = []; % Cannot select out of range of DEM, so remove those cols
    
    if wdir ~= 990 && wdir ~= 0 % wdir of 990 means wdir is variable, 0 means no wind
        orogliftarea = wspeed .* sind(slope(rows, cols)) .* cosd((wdir / 10) - aspect(rows, cols));
        oroglift = nanmax(orogliftarea(:));
    else
        oroglift = NaN;
    end
    
    wspeed_col(i) = wspeed;
    wdir_col(i) = wdir;
    oroglift_col(i) = oroglift;
    wstation_col(i) = cWS;
end

trackselection.oroglift = oroglift_col;
trackselection.wspeed = wspeed_col;
trackselection.wdir = wdir_col;
trackselection.wstation = wstation_col;

% Return annotated tracks
tracksOroLift = trackselection;

% Delete merged file
delete(merged_dem);

% % Now remove the initial tracks from the table
% tracks(tracks.longitude >= info.XMin & tracks.longitude <= info.XMax ...
%        & tracks.latitude >= info.YMin & tracks.latitude <= info.YMax, :) = [];
%    
% % And replace with the new ones containing data on orographic lift
% tracks = vertcat(tracks, trackselection);

end
