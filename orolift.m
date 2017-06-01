clear;
close all;
clc;

setenv('PATH', ['/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin', getenv('PATH')])
setenv('DYLD_LIBRARY_PATH',['/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin' getenv('DYLD_LIBRARY_PATH')])


profile on;

load('proj_settings_2015.mat');
load('proj_tracks_2015.mat');
load('proj_dems_2015.mat');
load('proj_wind_2015.mat');

% Prepare tracks to store orographic lift data
i = size(tracks, 1);
tracks.oroglift_max = zeros(i,1);
tracks.oroglift_mean = zeros(i,1);
tracks.oroglift_min = zeros(i,1);
tracks.dem_altitude = zeros(i,1);
tracks.wspeed = zeros(i,1);
tracks.wdir = zeros(i,1);
tracks.wstation = zeros(i,1);

n = size(dems, 1);

for i = 1:n
    tic
    deminfo = dems(i,:);
    
    % Select tracks within dem
    trackselection = tracks(tracks.longitude >= deminfo.XMin & tracks.longitude <= deminfo.XMax ...
                            & tracks.latitude >= deminfo.YMin & tracks.latitude <= deminfo.YMax, :);
    
    % Calculate orographic lift
    fprintf('%d/%d Starting orographic lift calculations on tile %s for %d records \n', i, n, string(deminfo.filenames), size(trackselection, 1));
    tracksOroLift = calcOrographicLift('data/dem/', deminfo.filenames, ...
                                       trackselection, stations, wind, ...
                                       2.5, 4, -1e10);
    
    % Remove old tracks without orographic lift data       
    tracks(tracks.longitude >= deminfo.XMin & tracks.longitude <= deminfo.XMax ...
           & tracks.latitude >= deminfo.YMin & tracks.latitude <= deminfo.YMax, :) = [];
       
    fprintf('%d/%d Calculated orographic lift for %d records on tile %s \n', i, n, size(tracksOroLift, 1), string(deminfo.filenames));
       
    % Now add tracks with orographic lift data
    tracks = vertcat(tracks, tracksOroLift);
    toc
end

save(['proj_tracks_oroglift_', daterange, '.mat'], 'tracks');

profile viewer;

%% Merge orographic lift data with classifications
csvfiles = dir(fullfile('data/classified/', '*.csv'));

classifications = [];
for i = 1:size(csvfiles,1)
    csv = readtable([csvfiles(i).folder '/' csvfiles(i).name]);
    classifications = vertcat(classifications, csv);
end

% Fix the formatting of some datetime elements
classifications.date_time = replace(classifications.date_time, 'T', ' ');
classifications.date_time = replace(classifications.date_time, '00Z', '');

% Merge with tracks to creat the final combinations of tracks information
% and behaviour classifications
classified_tracks = innerjoin(tracks, classifications);

save(['proj_tracks_classified_', daterange, '.mat'], 'classified_tracks');

%% Prepare data for data analysis
%  In our data analysis, we are only focussed on the flight strategies of
%  birds. We can therefore make a selection from the classified tracks on
%  rows that contain data on flying birds. Flight behaviour has been
%  classified with the following labels:
%
%  1. Flap
%  2. ExFlap
%  3. Soar
%  4. Boat
%  5. Float
%  6. SitStand
%  7. TerLoco
%  8. Other
%  9. Manouvre
%  10. Pecking
%
%  Classes 1, 2, 3 and 9 are flight-related, so these are the classes we
%  select for the data analysis.

flight_tracks = classified_tracks(classified_tracks.class_id == 1 | ...
                                  classified_tracks.class_id == 2 | ...
                                  classified_tracks.class_id == 3 | ...
                                  classified_tracks.class_id == 9, :);

% Save the final results in MATLAB and CSV file
save(['proj_flight_tracks_', daterange, '.mat'], 'flight_tracks');
writetable(flight_tracks, ['flight_tracks_', daterange, '.csv']);
