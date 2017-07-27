close all;
clc;

% ---------------------------------------------------------------------
% Research variables
% ---------------------------------------------------------------------

% Trackers to include in all queries
devices = csvread('data/devices.csv', 1, 0);
trackers_hg = sort(devices(devices(:,2) == 850));
trackers_lbbg = sort(devices(devices(:,2) == 15));
        
% Timeframe
startdatetime  = '2015-05-01 00:00:00';
enddatetime    = '2015-07-15 23:59:59';
datetimeformat = 'yyyy-mm-dd HH:MM:SS';

% Filenaming conventions
% Files are named with their corresponding year or year-range, finer
% control of naming could be set below
startyear = datestr(datenum(startdatetime, datetimeformat), 'yyyy');
endyear = datestr(datenum(enddatetime, datetimeformat), 'yyyy');
if strcmp(startyear, endyear)
    daterange = startyear;
else
    daterange = [startyear, '_', endyear];
end

% We only want high-resolution data, so we start by downloading
% accelerometer data, and subsequently download the corresponding tracker
% data only for those trackers which have provided a sufficiently large
% accelerometer dataset. Because the trackers for LBBG and HG differ, we
% have different minimum thresholds of datapoints per day (every datapoint
% corresponds with 20-40 accelerometer measurements). The choice for a
% minimum number of observations per day is somewhat arbitrary. In the case
% of LBBG, there are so few records that we use all data and thus the
% minimum # of observations is set to 0

hg_min_obs = 2880; % 'once every 30 seconds' (not really)
lbbg_min_obs = 0; % Too low data resolution, so we accept everything

%% Fetch Accelerometer data
highres_hg_devices = prepBirdAccelerometerData(db_user, db_pass, ...
    trackers_hg, startdatetime, enddatetime, hg_min_obs, ...
    'data/ResearchAreaNH.shp', ['data/unclassified/', daterange]);

highres_lbbg_devices = prepBirdAccelerometerData(db_user, db_pass, ...
    trackers_lbbg, startdatetime, enddatetime, lbbg_min_obs, ...
    'data/ResearchAreaNH.shp', ['data/unclassified/', daterange]);

save(['temp/proj_settings_', daterange, '.mat'], 'db_pass', 'db_user', 'devices', ...
     'enddatetime', 'hg_min_obs', 'highres_hg_devices', ...
     'highres_lbbg_devices', 'lbbg_min_obs', 'startdatetime', ...
     'trackers_hg', 'trackers_lbbg', 'daterange');
 
%% Prepare accelerometer data files for classification
%  In order to classify based on the accelerometer, do the following: 
%  1. Run the code in this section to generate a paths.txt file containing
%     paths to all the .mat files with accelerometer data. Subsequently
%     copy all the .mat files from the data/unclassified/[year] over to the
%     classifier.
%  2. Copy the paths from paths.txt over to the settings.properties file of
%     classifier. All lines are commented out by default. The
%     fastest/easiest way to do the classification is by uncommenting these
%     lines one by one and running the classifier on each one of these
%     lines. Every time a classifications.csv file is created, you can
%     rename these to the respective device ID. So after classification,
%     you end up with files such as 317.csv, 6202.csv, etc.
%
%  @Note: Sometimes (very rarely) accelerometer data contains NaNs in the
%  tspeed columns (and possibly other columns as well), these cannot be
%  classified by the classifier and will throw an error. Possibly this can
%  be fixed, but for now you will have to identify this file manually and
%  change things manually. This is a bit burdensome, but can be done by
%  feeding successively halves of the dataset to the classifier and
%  checking when it throws an error. If it does, feed half of that data
%  again until you find the culprit, etc. In this file all you need to do
%  is remove the NaN records.

matfiles = dir(fullfile(['data/unclassified/', daterange, '/'], '*.mat'));

filenames = {matfiles(:).name}';

device_ids = strtok(filenames, 'd'); % d is the delimiter between ID and yyyymmdd
unique_ids = unique(device_ids);
n = numel(unique_ids);

% Now build these into classifier-compatible strings
lines = {};

for i = 1:n
    % Select all filenames with a similar ID prefix
    current_id = unique_ids(i);
    selector = strcmp(device_ids, current_id);
    fileselection = filenames(selector);
    
    % Now transform that in a string
    paths = join(fileselection, ',');
    space = {' = '};
    lines{i,1} = strcat('# unannotated_measurement_source_paths', space, paths);
end

fileID = fopen(['data/unclassified/', daterange, '/paths.txt'], 'w');
fprintf(fileID, '%s \n', lines{:});
fclose(fileID);

 
%% Fetch normal tracking data
%  Now that we now which birds have datapoints with a sufficient resolution
%  (excluding lbbg) and data within the given timeframe, we can download
%  the tracker data from these moments
trackers = [highres_hg_devices, highres_lbbg_devices];
tracks = prepBirdData(db_user, db_pass, trackers, startdatetime, enddatetime, 'data/ResearchAreaNH.shp');
save(['temp/proj_tracks_', daterange, '.mat'], 'tracks');

%% Fetch weather (wind) data
%  In order to calculate orographic lift, we need wind data for all
%  locations within the shapefile.
startdatenum = datenum(startdatetime, datetimeformat);
enddatenum   = datenum(enddatetime, datetimeformat);
conversionformat = 'yyyymmdd';
startdate = str2double(datestr(startdatenum, conversionformat));
enddate   = str2double(datestr(enddatenum, conversionformat));

[wind, stations] = prepWindData(['data/weather/', daterange, '/'], 'data/ResearchAreaNH.shp', startdate, enddate);
wind = table(wind(:,1), wind(:,2), wind(:,3), wind(:,4), wind(:,5), wind(:,6), wind(:,7), ...
             'VariableNames', {'stationID', 'date', 'hour', 'wdir', 'wspeed_hr', 'wspeed_10min', 'wspeed_peak'});
save(['temp/proj_wind_', daterange, '.mat'], 'wind', 'stations');

%% Load DEM info from files
%  For the sake of simplicity, we are not going to seperately download the
%  DEM files again. Code is properly documented in the prepDEM file. Now we
%  are only making a table containing the metadata of the DEM tiles.
dems = loadDEMinfo('data/dem/', '.wgs84.tif', -1e10);
save(['temp/proj_dems_', daterange, '.mat'], 'dems');
