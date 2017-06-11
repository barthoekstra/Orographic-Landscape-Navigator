clear;
close all;
clc;

addpath(genpath('functions/'));

tracks = load('data/flight_tracks_2015_2016.mat');
tracks = tracks.flight_tracks;

%% Add altitude above DEM level column to dataset
% Since we assume orographic lift soaring happens close to the surface,
% this is a relevant column to use to select the right datapoints.
tracks.altitude_adl = tracks.altitude - tracks.dem_alt_mean;

%% Select soaring datapoints
% In order to determine the characteristics of orographic lift soaring, we
% need to make a selection of high-resolution tracks of soaring flight,
% because those offer the most accurate data. We thus:
% 1. Select soaring tracks with consecutively more than 5 datapoints which
%    are classified as soaring by the classifier.
% 2. Select only datapoints where the altitude above the mean DEM altitude
%    in meters is higher than 0 and lower than 25, because we assume
%    orographic lift soaring happens fairly close to the surface
consecutive_min = 5;
adl_min = 0; % adl = above dem level
adl_max = 25;

soaring = tracks(tracks.class_id == 3 & tracks.altitude_adl > adl_min & ...
                 tracks.altitude_adl <= adl_max, :);
             
soaring_all = tracks(tracks.class_id == 3, :);
             
%% Determine consecutive track lengths
% In order to find high-resolution and highly accurate tracks, we determine
% the length of consecutively increasing ID integers. So [1, 2, 3, 4, 6]
% result in [4, 1]. The following code calculates these values, which we
% can subsequently use to select long tracks of consecutive numbers.
IDs = soaring.ID;
consecutive_ints = diff(find(diff([nan; IDs(:); nan]) ~= 1));
consecutive_expanded = rude(consecutive_ints, consecutive_ints);
soaring.consecutive = consecutive_expanded';
soaring = [soaring(:,1), soaring(:,40), soaring(:,2:39)];

%% Selection of tracks for analysis
% Now that we have all soaring datapoints between 0 and 25 meters above the
% DEM level and their number of consecutive measurements, we can select
% only the tracks that have a lengt of 5 consecutive datapoints.
%
% Note: Following this methodology it is theoretically possible two tracks
% of two different birds are taken as one, so build in functionality that
% checks for this
hr_soaring = soaring(soaring.consecutive >= consecutive_min, :);

%% Calculate summary statistics for these soaring tracks
% We want to gather the following for every soaring track
% 1. altitude_adl (max, min, mean and std)
% 2. oroglift (max, min, mean and std) (based on oroglift_max values)
% 3. start and stop IDs and track_length
% 4. tracker device_info_serial
% 5. begin time and end time
% 6. project

% Determine # of tracks
% differences = diff(hr_soaring.consecutive);
% nr_tracks = sum(differences ~= 0);
% col = NaN(nr_tracks, 1);

% Preallocate arrays
alt_adl_max = []; alt_adl_min = []; alt_adl_mean = []; alt_adl_std = [];
oroglift_max = []; oroglift_min = []; oroglift_mean = []; oroglift_std = [];
start_id = []; stop_id = []; track_length = [];
device_info_serial = [];
begin_datetime = {}; end_datetime = {};
project = {};
sel_len = [];
wdir = []; wspeed = [];

n = size(hr_soaring, 1);

k = 1;
i = 1;
while i <= n
    selection_length = hr_soaring(i,:).consecutive;
    
    % Now select the track, but remove end datapoints, because these
    % values, may not accurately represent the conditions for soaring,
    % since the bird decides to stop soaring after these datapoints
    track = hr_soaring(i:(i + selection_length - 1), :);
    
    % Store summary statistics
    alt_adl_max(k) = max(track.dem_alt_mean);
    alt_adl_min(k) = min(track.dem_alt_mean);
    alt_adl_mean(k) = mean(track.dem_alt_mean);
    alt_adl_std(k) = std(track.dem_alt_mean);
    
    oroglift_max(k) = max(track.oroglift_max);
    oroglift_min(k) = min(track.oroglift_max);
    oroglift_mean(k) = mean(track.oroglift_max);
    oroglift_std(k) = std(track.oroglift_max);
    
    % Track information
    start_id(k) = track(1,:).ID;
    stop_id(k) = track(end,:).ID;
    track_length(k) = stop_id(k) - start_id(k);
    device_info_serial(k) = unique(track.device_info_serial);
    begin_datetime(k) = track(1,:).date_time;
    end_datetime(k) = track(end,:).date_time;
    project(k) = unique(track.project);
    wdir(k) = mean(track.wdir);
    wspeed(k) = mean(track.wspeed);
    
    % Other
    sel_len(k) = selection_length;
    
    k = k + 1;
    i = i + selection_length;
end

varnames = {'device_info_serial', 'start_id', 'stop_id', 'begin_datetime', ...
            'end_datetime', 'track_length', 'project', 'oroglift_max', ...
            'oroglift_mean', 'oroglift_min', 'oroglift_std', 'alt_adl_max', ...
            'alt_adl_mean', 'alt_adl_min', 'alt_adl_std', 'wdir', 'wspeed'};

summary = table(device_info_serial', start_id', stop_id', begin_datetime', ...
                        end_datetime', track_length', project', oroglift_max', ...
                        oroglift_mean', oroglift_min', oroglift_std', ...
                        alt_adl_max', alt_adl_mean', alt_adl_min', alt_adl_std', ...
                        wdir', wspeed', ...
                        'VariableNames', varnames);
                    
% We can also see there are some tracks occurring that are over the sea,
% since alt_adl_max equals zero. Orographic lift cannot be generated in
% these locations, so we remove these rows from the dataset as well.
summary = summary(summary.alt_adl_max ~= 0, :);
                    
%% Prepare randomized tracks to be included in the analysis
% We have also made random tracks, where all values are equal, but lat and
% lon values have been randomized, to simulate the effect of birds moving
% around without having or using any knowledge about the characteristics of
% the area. In other words: they 'fly' around as if they do not know about
% the existence of thermal lift, orographic lift or strategies that make
% flying easier/less costly in terms of energy expenditure.
random_tracks = load('data/random_tracks_2015_2016.mat');
random_tracks = random_tracks.random_tracks;

% Similar to the real tracks, we remove the tracks over the sea
random_tracks = random_tracks(random_tracks.dem_alt_max ~= 0, :);

%% Make subselections for different species
soaring_hg = soaring(strcmp(soaring.project, 'HG_TEXEL') == 1, :);
soaring_lbbg = soaring(strcmp(soaring.project, 'LBBG_TEXEL') == 1, :);
hr_soaring_hg = hr_soaring(strcmp(hr_soaring.project, 'HG_TEXEL') == 1, :);
hr_soaring_lbbg = hr_soaring(strcmp(hr_soaring.project, 'LBBG_TEXEL') == 1, :);
%% Compare soaring with random flight patterns
% Both species
pd_soaring = fitdist(soaring.oroglift_max, 'Kernel');
pd_soaring_consec = fitdist(hr_soaring.oroglift_max, 'Kernel');
pd_random = fitdist(random_tracks.oroglift_max, 'Kernel');

% Herring gull
pd_soaring_hg = fitdist(soaring_hg.oroglift_max, 'Kernel');
pd_soaring_hg_consec = fitdist(hr_soaring_hg.oroglift_max, 'Kernel');
pd_soaring_lbbg = fitdist(soaring_lbbg.oroglift_max, 'Kernel');
pd_soaring_lbbg_consec = fitdist(hr_soaring_lbbg.oroglift_max, 'Kernel');

x = -2:.1:12;
yrandom = pdf(pd_random, x);

y1 = pdf(pd_soaring, x); % Both species of gulls combined
y2 = pdf(pd_soaring_consec, x);
y1hg = pdf(pd_soaring_hg, x); % Herring Gull
y2hg = pdf(pd_soaring_hg_consec, x);
y1lbbg = pdf(pd_soaring_lbbg, x); % Lesser Black-backed Gull
y2lbbg = pdf(pd_soaring_lbbg_consec, x);

% Where are the intersections? In terms of orographic lift:
[x_soaring_all, ~] = intersections(x, yrandom, x, y1, 1)
[x_soaring_consec, ~] = intersections(x, yrandom, x, y1, 1)
[x_soaring_hg, ~] = intersections(x, yrandom, x, y1hg, 1)
[x_soaring_hg_consec, ~] = intersections(x, yrandom, x, y2hg, 1)
[x_soaring_lbbg, ~] = intersections(x, yrandom, x, y1lbbg, 1)
[x_soaring_lbbg_consec, ~] = intersections(x, yrandom, x, y2lbbg, 1)

figure;
plot(x, y2hg, 'DisplayName', 'Long soaring (Herring Gull)', 'LineWidth', 2);
hold on;
plot(x, y2lbbg, 'DisplayName', 'Long soaring (Lesser Black-backed Gull)', 'LineWidth', 2);
plot(x, yrandom, 'DisplayName', 'Background lift', 'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
legend('show');
title('Orographic lift rates: Herring vs Lesser Black-backed Gull');
ylabel('Probability density');
xlabel('Orographic lift [m/s]');

figure;
plot(x, y1hg, 'DisplayName', 'Soaring Herring Gulls', 'LineWidth', 2);
hold on;
plot(x, y1lbbg, 'DisplayName', 'Soaring Lesser Black-backed Gulls', 'LineWidth', 2);
plot(x, yrandom, 'DisplayName', 'Background lift', 'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
legend('show');
title('Orographic lift rates: Herring vs Lesser Black-backed Gull');
ylabel('Probability density');
xlabel('Orographic lift [m/s]');

%%
oroglift_lower_bound = 0.8141;
oroglift_upper_bound = 4.2078;
orog_soaring = tracks(tracks.oroglift_max >= oroglift_lower_bound & ...
                      tracks.oroglift_max <= oroglift_upper_bound & ...
                      tracks.altitude_adl < 100 & ...
                      tracks.altitude_adl > 0 & ...
                      tracks.class_id == 3, :);
%% Influence of wind direction and landscape aspect
pd_landscape_soaring = fitdist(hr_soaring.landscape, 'Kernel');
pd_landscape_flapping = fitdist(flapping.landscape, 'Kernel');
x = -10:1:360;
y1 = pdf(pd_landscape_soaring, x);
y2 = pdf(pd_landscape_flapping, x);
figure;
plot(x, y1);
hold on;
plot(x, y2);