clear;
close all;
clc;

setenv('PATH', ['/usr/local/bin:/Users/barthoekstra/anaconda/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/Library/TeX/texbin:', getenv('PATH')])

profile on;

load('proj_settings.mat');
load('proj_tracks.mat');
load('proj_dems.mat');
load('proj_wind.mat');

% Prepare tracks to store orographic lift data
i = size(tracks, 1);
tracks.oroglift = zeros(i,1);
tracks.wspeed = zeros(i,1);
tracks.wdir = zeros(i,1);
tracks.wstation = zeros(i,1);

n = size(dems, 1);

for i = 1:n
    deminfo = dems(i,:);
    
    % Select tracks within dem
    trackselection = tracks(tracks.longitude >= deminfo.XMin & tracks.longitude <= deminfo.XMax ...
                            & tracks.latitude >= deminfo.YMin & tracks.latitude <= deminfo.YMax, :);
    
    % Calculate orographic lift 
    tracksOroLift = calcOrographicLift('data/dem/', deminfo.filenames, ...
                                       trackselection, stations, wind, ...
                                       2.5, 4, -1e10);
    
    % Remove old tracks without orographic lift data       
    tracks(tracks.longitude >= deminfo.XMin & tracks.longitude <= deminfo.XMax ...
           & tracks.latitude >= deminfo.YMin & tracks.latitude <= deminfo.YMax, :) = [];
       
    fprintf('%d/%d Calculated orographic lift for %d records on tile %s \n', i, n, size(tracksOroLift, 1), string(deminfo.filenames));
       
    % Now add tracks with orographic lift data
    tracks = vertcat(tracks, tracksOroLift);
end

profile viewer;

% Speed:
% Without parfor in calcOrographicLift loop: 978s (16 min)
% With parfor in calcOrographicLift loop: 