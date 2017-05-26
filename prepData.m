close all;
clc;

% Time to define some research variables. 

% Trackers to include in all queries
devices = csvread('data/devices.csv', 1, 0);
trackers_hg = sort(devices(devices(:,2) == 850));
trackers_lbbg = sort(devices(devices(:,2) == 15));
        
% Timeframe
startdatetime  = '2015-05-01 00:00:00';
enddatetime    = '2015-07-15 23:59:59';
datetimeformat = 'yyyy-mm-dd HH:MM:SS';

% We only want high-resolution data, so we start by downloading
% accelerometer data, and subsequently download the corresponding tracker
% data only for those trackers which have provided a sufficiently large
% accelerometer dataset. Because the trackers for LBBG and HG differ, we
% have different minimum thresholds of datapoints per day (every datapoint
% corresponds with 20-40 accelerometer measurements). The choice for a
% minimum number of observations per day is somewhat arbitrary. In the case
% of LBBG, there are so few records that we use all data and thus the
% minimum # of observations is set to 0
% @NOTE: In further calculations, effects have to be corrected for the
% number of observations in a day (!)
hg_min_obs = 2880; % 'once every 30 seconds' (not really)
lbbg_min_obs = 0; % Usage of 180 results in data for 606, 754, 805, 806 and 871

%% Fetch Accelerometer data
highres_hg_devices = prepBirdAccelerometerData(db_user, db_pass, ...
    trackers_hg, startdatetime, enddatetime, hg_min_obs, ...
    'data/ResearchAreaNH.shp', 'data/unclassified');

highres_lbbg_devices = prepBirdAccelerometerData(db_user, db_pass, ...
    trackers_lbbg, startdatetime, enddatetime, lbbg_min_obs, ...
    'data/ResearchAreaNH.shp', 'data/unclassified');

save('proj_settings.mat', 'db_pass', 'db_user', 'devices', ...
     'enddatetime', 'hg_min_obs', 'highres_hg_devices', ...
     'highres_lbbg_devices', 'lbbg_min_obs', 'startdatetime', ...
     'trackers_hg', 'trackers_lbbg');
%% Fetch normal tracking data
%  Now that we now which birds have datapoints with a sufficient resolution
%  (excluding lbbg) and data within the given timeframe, we can download
%  the tracker data from these moments
trackers = [highres_hg_devices, highres_lbbg_devices];
tracks = prepBirdData(db_user, db_pass, trackers, startdatetime, enddatetime, 'data/ResearchAreaNH.shp');
save('proj_tracks.mat', 'tracks');

%% Fetch weather (wind) data
%  In order to calculate orographic lift, we need wind data for all
%  locations within the shapefile.
startdatenum = datenum(startdatetime, datetimeformat);
enddatenum   = datenum(enddatetime, datetimeformat);
conversionformat = 'yyyymmdd';
startdate = str2double(datestr(startdatenum, conversionformat));
enddate   = str2double(datestr(enddatenum, conversionformat));

[wind, stations] = prepWeatherData('data/weather/', 'data/ResearchAreaNH.shp', startdate, enddate);
save('proj_wind.mat', 'wind', 'stations');

%% Load DEM info from files
%  For the sake of simplicity, we are not going to seperately download the
%  DEM files again. Code is properly documented in the prepDEM file, but
%  gdal and ogr2ogr commands may have to be called from Python, for which a
%  script is still to be written
deminfo = loadDEMinfo('data/dem/', '.wgs84.tif', -1e10);
save('proj_deminfo.mat', 'deminfo');