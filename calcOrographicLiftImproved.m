%function [ output_args ] = calcOrographicLift( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%   setenv('PATH', ['/usr/local/bin:/Users/barthoekstra/anaconda/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/Library/TeX/texbin:', getenv('PATH')])

% Some global variables
profile on;
nanthreshold = -1e10;
cellsize = 2.5; % m (meter)
stepsize = 4; % 10m (meter)
storpath = 'data/dem/';
filepostfix = '.wgs84.tif';

% Load DEM and change values under threshold value to be NaN, so they do
% not affect our calculations
%filename = 'r14bn2.wgs84.tif';
filename = 'r19az1.wgs84.tif';

% ---------------------------------------------------------------------
% Prepare tracks table to contain relevant orographic lift information
% ---------------------------------------------------------------------
num_tracks = size(tracks, 1);
tracks.oroglift = zeros(num_tracks, 1);
tracks.wspeed   = tracks.oroglift;
tracks.wdir     = tracks.oroglift;
tracks.wstation = tracks.oroglift;

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
ref_tile_name = replace(filename,'.wgs84.tif','');
ref_tile_name = ref_tile_name(2:end);

ref_tile = tiles(strcmp(ref_tile_name, tiles.UNIT),:);

merge_units = [];

for i = 1:size(tiles,1)
   comp_tile = tiles(i,:);
   
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
files = join(strcat(storpath, 'r', merge_units, filepostfix), ' ');
output = strcat(storpath, 'r', ref_tile_name, '.merged.wgs84.tif');
cmd = sprintf('python gisactions.py merge %s %s', output, files);
          
[status, cmd] = system(cmd);
if status == 1
    error(cmd)
end

%%
% ---------------------------------------------------------------------
% Select all datapoints within boundary of the original (non-merged) 
% DEM
% ---------------------------------------------------------------------
original_dem = [storpath, filename];
merged_dem = output;
info = deminfo(strcmp(deminfo.filenames, filename),:);

trackselection = tracks(tracks.longitude >= info.XMin & tracks.longitude <= info.XMax ...
                        & tracks.latitude >= info.YMin & tracks.latitude <= info.YMax, :);

%%
% ---------------------------------------------------------------------
% Calculate orographic lift
% ---------------------------------------------------------------------
dem = geotiffread(output);
dem(dem <= nanthreshold) = 0;

[r, c] = size(dem);

% dzdy
% dzdy = zeros(r, c);
% dzdy(1,:) = (dem(2,:) - dem(1,:)) ./ (cellsize);
% dzdy(r,:) = (dem(r,:) - dem(r-1,:)) ./ (cellsize);
% dzdy(2:r-1,:) = (dem(3:r,:) - dem(1:r-2,:)) ./ (cellsize * 2);
dzdy = zeros(r, c);
dzdy(1,:) = (dem(1,:) - dem(2,:)) ./ (cellsize);
dzdy(r,:) = (dem(r-1,:) - dem(r,:)) ./ (cellsize);
dzdy(2:r-1,:) = (dem(1:r-2,:) - dem(3:r,:)) ./ (cellsize * 2);

% dzdx
% dzdx = zeros(r, c);
% dzdx(:,1) = (dem(:,2) - dem(:,1)) ./ (cellsize);
% dzdx(:,c) = (dem(:,c) - dem(:,c-1)) ./ (cellsize);
% dzdx(:,2:c-1) = (dem(:,3:c) - dem(:,1:c-2)) ./ (cellsize * 2);
dzdx = zeros(r, c);
dzdx(:,1) = (dem(:,1) - dem(:,2)) ./ (cellsize);
dzdx(:,c) = (dem(:,c-1) - dem(:,c)) ./ (cellsize);
dzdx(:,2:c-1) = (dem(:,1:c-2) - dem(:,3:c)) ./ (cellsize * 2);

% Slope
slope = atand(sqrt(dzdx.^2 + dzdy.^2));

% Aspect
% aspect = 90 + atan2d(dzdy, dzdy);
% aspect = 90 + atan2d(dzdy, dzdx);
aspect = 90 + atan2d(dzdy, dzdx);
aspect = aspect + (aspect < 0) * 360; % transform to [0 360]

% Cleanup
%aspect(isnan(aspect)) = 0; % @NOTE: Is it right to set this to 0? Probably not, but alternative?
slope(isnan(slope)) = 0;

% Calculate Ca (Updraft coefficient) for every 10 degrees change of wind
% direction. We don't need to calculate Ca in finer steps, because KNMI
% only reports the values in steps of 10 degrees.
Ca = zeros(size(dem, 1), size(dem, 2), 36);
for i = 10:10:360
    k = i / 10;
    Ca(:,:,k) = sind(slope) .* cosd(i - aspect);
end

%%
% Now it is - finally - time to really calculate the orographic lift
spatialinfo = geotiffinfo(merged_dem);
n = size(trackselection,1);

% Make some variables to store column data in, so we can merge them with
% the table later on
wspeed_col   = zeros(n,1);
wdir_col     = wspeed_col;
oroglift_col = wspeed_col;
wstation_col = wspeed_col;

for i = 1:n
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
    rows = r-10:1:r+10;
    cols = c-10:1:c+10;
    if wdir ~= 990
        orogliftarea = wspeed .* Ca(rows, cols, wdir / 10);
        oroglift = max(orogliftarea(:));
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

tracks(tracks.longitude >= info.XMin & tracks.longitude <= info.XMax ...
       & tracks.latitude >= info.YMin & tracks.latitude <= info.YMax, :) = ...
       trackselection;
%end

%%
% %%
% dem = geotiffread([storpath, filename]);
% dem(dem <= nanthreshold) = NaN;
% info = deminfo(strcmp(deminfo.filenames, filename),:)
% spatialref = geotiffinfo(['data/dem/', filename]);
% boundingbox = spatialref.BoundingBox;
% 
% % Now select all the datapoints within the bounds of this DEM
% trackselection = tracks(tracks.longitude >= info.XMin & tracks.longitude <= info.XMax ...
%                         & tracks.latitude >= info.YMin & tracks.latitude <= info.YMax, :);
% 
% % Now we calculate the slope and aspect of the entire DEM
% 
% % dzdx - Change in Easterly direction
% dzdx = (dem(:,1:end-2) - dem(:,3:end)) ./ (cellsize * 2);
% dzdxleft = (dem(:,1) - dem(:,2)) ./ (cellsize);
% dzdxright = (dem(:,end-1) - dem(:,end)) ./ (cellsize);
% dzdx = [dzdxleft, dzdx, dzdxright];
% 
% % dzdy - Change in Northerly direction
% dzdy = (dem(1:end-2,:) - dem(3:end,:)) ./ (cellsize * 2);
% dzdytop = (dem(1,:) - dem(2,:)) ./ (cellsize);
% dzdybottom = (dem(end-1,:) - dem(end,:)) ./ (cellsize);
% dzdy = [dzdytop; dzdy; dzdybottom];
% 
% % Slope
% slope = atand(sqrt(dzdx.^2 + dzdy.^2));
% 
% % Aspect
% aspect = 90 + atan2d(dzdy, dzdx);
% aspect = aspect + (aspect < 0) * 360; % transform to [0 360]
% 
% % Cleanup
% aspect(isnan(aspect)) = 0; % @NOTE: Is it right to set this to 0? Probably not, but alternative?
% slope(isnan(slope)) = 0;
% 
% % Now calculate the updraft coefficient for every possible wind direction
% % excluding no wind (wdir = 0) and variable wind (wdir = 990). Since wind
% % direction is only logged in steps of 10 deg, we only have to calculate
% % the updraft coefficient 36 times, once for every step of 10 degrees in
% % 360 degrees
% for i = 10:10:360
%     k = i / 10;
%     Ca(:,:,k) = sind(slope) .* cosd(i - aspect);
% end
% 
% figure;
% imagesc(linspace(boundingbox(1), boundingbox(2), spatialref.Width), ...
%         linspace(boundingbox(3), boundingbox(4), spatialref.Height), ...
%         flipud(dem));
% set(gca,'YDir','normal');
% drawnow;
% hold on;
%%
% Now it is - finally - time to really calculate the orographic lift
% i = 1;
% for i = 1:size(trackselection,1)
%     datetime = datenum(trackselection(i,:).date_time, 'yyyy-mm-dd HH:MM:SS');
%     date = str2num(datestr(datetime, 'yyyymmdd'));
%     tracksel = trackselection(i,:);
%     lat = tracksel.latitude;
%     lon = tracksel.longitude;
%     windselection = wind(wind.stationID == 209 & wind.date == date ...
%                          & wind.hour == trackselection(i,:).hour,:);
%     wspeed = windselection.wspeed_hr / 10;
%     wdir   = windselection.wdir / 10;
%     
%     [r, c] = setpostn(dem, spatial.SpatialRef, lat, lon);
%     oroglift = wspeed * Ca(r, c, wdir);
%     
%     plot(lon, lat, 'r', 'Marker', 'o', 'MarkerSize', 10);
%     xlim([boundingbox(1), boundingbox(2)]); 
%     ylim([boundingbox(3), boundingbox(4)]);
% end
% 
% %end
% 

%% Orog Lift works odd, let's check
% Catest = Ca(2800:4300,2200:3000,:);
% 
% figure;
% subplot(2,2,1);
% imagesc(Catest(:,:,12)); colorbar; colormap(bluewhitered)
% title('wdir: 120');
% subplot(2,2,2);
% imagesc(Catest(:,:,18)); colorbar; colormap(bluewhitered)
% title('wdir: 180');
% subplot(2,2,3);
% imagesc(Catest(:,:,27)); colorbar; colormap(bluewhitered)
% title('wdir: 270');
% subplot(2,2,4);
% imagesc(Catest(:,:,36)); colorbar; colormap(bluewhitered)
% title('wdir: 360/0');

% %% Plot orog lift
% orig_info = geotiffinfo(original_dem);
% orig_dem = geotiffread(original_dem);
% orig_dem(orig_dem <= nanthreshold) = 0;
% 
% figure;
% %subplot(2,2,1);
% imagesc([info.XMin info.XMax], [info.YMin info.YMax], flipud(orig_dem));
% colormap(bluewhitered);
% set(gca, 'YDir', 'normal');
% hold on;
% %selection = trackselection(trackselection.wdir >= 200 & trackselection.wdir <= 220,:);
% %selection = trackselection(trackselection.wdir >= 10 & trackselection.wdir <= 30,:);
% selection = trackselection(trackselection.wdir >= 110 & trackselection.wdir <= 130,:);
% plot(selection.longitude, selection.latitude, 'ok', 'MarkerSize', 5)
% xlim([info.XMin info.XMax]); ylim([info.YMin info.YMax]);
% %%
% hold on;
% scatter(trackselection.longitude, trackselection.latitude, 15, trackselection.oroglift);
% colormap('parula');


profile viewer;


