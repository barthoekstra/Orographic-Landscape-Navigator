clear;
close all;
clc;
%%
addpath(genpath('functions/'));

tracks = load('data/flight_tracks_2015_2016.mat');
tracks = tracks.flight_tracks;

% Also check for night results? Set to 1
nightcheck = 1;
nighthours = [21 23 0 4]; % hour 21-23 of the night's start day, hour 0-4 of the night's end day

% Make sure the table is sorted
tracks = sortrows(tracks,{'year', 'ID'},{'ascend','ascend'});

%% Add altitude above DEM level column to dataset
% Since we assume orographic lift soaring happens close to the surface,
% this is a relevant column to use to select the right datapoints.
tracks.altitude_adl = tracks.altitude - tracks.dem_alt_mean;

%% Only select tracks with consecutively the same flight strategy class
IDs = tracks.ID;
consecutive_ints = diff(find(diff([nan; IDs(:); nan]) ~= 1));
consecutive_expanded = rude(consecutive_ints, consecutive_ints);
tracks.consecutive = consecutive_expanded';
%tracks = tracks(tracks.consecutive >= 5, :);

%% Summarize tracks
devices = unique(tracks.device_info_serial);
n_devices = size(devices, 1);
years = unique(tracks.year);
n_years = size(years, 1);

% Orographic lift parameters
max_meters = 30;
oroglift_min_all = 0.8141; % minimum value hg + lbbg combined
oroglift_max_all = 4.2078; % maximum value hg + lbbg combined
oroglift_min_hg = 0.7695; % minimum value hg
oroglift_max_hg = 4.1916; % maximum value hg
oroglift_min_lbbg = 0.8505; % minimum value lbbg
oroglift_max_lbbg = 4.3115; % maximum value lbbg

% Set the values used for the analysis below here:
oroglift_min = oroglift_min_all;
oroglift_max = oroglift_max_all;

% @ REMOVE THIS
%oroglift_min = 0.4

% Prepare columns
n_obs_year = [];
n_flight_obs_year = [];
n_flight_night_obs_year = [];
n_soaring_obs_year = [];
n_soaring_night_obs_year = [];
n_orog_soaring_all_year = [];
n_orog_soaring_night_all_year = [];
n_orog_soaring_obs_year = zeros(1, max_meters);
n_orog_soaring_night_obs_year = zeros(1, max_meters);
device = [];
year = [];
project = [];

for i = 1:n_devices
    for j = 1:n_years
        flight_sel = tracks(tracks.device_info_serial == devices(i) & tracks.year == years(j), :);
        
        if size(flight_sel, 1) == 0
            continue;
        end
        
        year = [year; years(j)];
        device = [device; devices(i)];
        project = [project; string(flight_sel(1,:).project)];
        
        summarystats = grpstats(flight_sel, {'month', 'day'}, 'mean', 'DataVars', 'n_obs_day');
        
        n_obs_year = [n_obs_year; sum(summarystats.mean_n_obs_day)];
        n_flight_obs_year = [n_flight_obs_year; size(flight_sel, 1)];
        
        % Night flight
        if nightcheck == 1
            night_sel = flight_sel((flight_sel.hour >= nighthours(1) & ...
                                    flight_sel.hour <= nighthours(2)) | ...
                                    (flight_sel.hour >= nighthours(3) & ...
                                    flight_sel.hour <= nighthours(4)), :);
            n_flight_night_obs_year = [n_flight_night_obs_year; size(night_sel, 1)];
        end
        
        % Select soaring
        soaring_sel = flight_sel(flight_sel.class_id == 3, :);
        n_soaring_obs_year = [n_soaring_obs_year; size(soaring_sel, 1)];
        
        % Select night soaring
        if nightcheck == 1
            night_soaring_sel = flight_sel(flight_sel.class_id == 3 & ...
                                           ((flight_sel.hour >= nighthours(1) & ...
                                           flight_sel.hour <= nighthours(2)) | ...
                                           (flight_sel.hour >= nighthours(3) & ...
                                           flight_sel.hour <= nighthours(4))), :);
            n_soaring_night_obs_year = [n_soaring_night_obs_year; size(night_soaring_sel, 1)];
        end
        
        % Select orographic soaring without altitude
        % Criteria:
        % 1. Oroglift over 0.9m/s and below 4.2m/s
        orog_soaring_sel = soaring_sel(soaring_sel.oroglift_max >= oroglift_min & ...
                           soaring_sel.oroglift_max <= oroglift_max, :);
        n_orog_soaring_all_year = [n_orog_soaring_all_year; size(orog_soaring_sel, 1)];
        
        % Select nightly orographic soaring
        if nightcheck == 1
            orog_night_soaring_sel = flight_sel(flight_sel.class_id == 3 & ...
                                               flight_sel.oroglift_max >= oroglift_min & ...
                                               flight_sel.oroglift_max <= oroglift_max & ...
                                               ((flight_sel.hour >= nighthours(1) & ...
                                               flight_sel.hour <= nighthours(2)) | ...
                                               (flight_sel.hour >= nighthours(3) & ...
                                               flight_sel.hour <= nighthours(4))), :);
            n_orog_soaring_night_all_year = [n_orog_soaring_night_all_year; size(orog_night_soaring_sel, 1)];
        end
        
        % Select orographic soaring with altitude
        % Criteria:
        % 1. Oroglift >= oroglift_min and <= oroglift_max
        % 2. Altitude_adl <= 1:max_meters
        orog_soaring = [];
        for k = 1:max_meters
            orog_soaring_sel = soaring_sel(soaring_sel.oroglift_max >= oroglift_min & ...
                                       soaring_sel.oroglift_max <= oroglift_max & ...
                                       soaring_sel.altitude_adl <= k, :);
            orog_soaring(1,k) = size(orog_soaring_sel, 1);
        end
        n_orog_soaring_obs_year = [n_orog_soaring_obs_year; orog_soaring];
        
        if nightcheck == 1
            orog_soaring_night = [];
            for k = 1:max_meters
                orog_soaring_sel = flight_sel(flight_sel.class_id == 3 & ...
                                              flight_sel.oroglift_max >= oroglift_min & ...
                                              flight_sel.oroglift_max <= oroglift_max & ...
                                              flight_sel.altitude_adl <= k & ...
                                              ((flight_sel.hour >= nighthours(1) & ...
                                              flight_sel.hour <= nighthours(2)) | ...
                                              (flight_sel.hour >= nighthours(3) & ...
                                              flight_sel.hour <= nighthours(4))), :);
                orog_soaring_night(1,k) = size(orog_soaring_sel, 1);
            end
            n_orog_soaring_night_obs_year = [n_orog_soaring_night_obs_year; orog_soaring_night];
        end
        
    end
end

n_orog_soaring_obs_year = n_orog_soaring_obs_year(2:end, :);
n_orog_soaring_night_obs_year = n_orog_soaring_night_obs_year(2:end, :);

varnames = {'device', 'year', 'project', 'n_obs', 'n_flight', ...
            'n_soaring', 'n_oroglift_all', 'n_oroglift'};
        
quantified = table(device, year, project, n_obs_year, n_flight_obs_year, ...
                   n_soaring_obs_year, n_orog_soaring_all_year, n_orog_soaring_obs_year, ...
                   'VariableNames', varnames);
               
quantified_night = table(device, year, project, n_obs_year, n_flight_night_obs_year, ...
                        n_soaring_night_obs_year, n_orog_soaring_night_all_year, n_orog_soaring_night_obs_year, ...
                        'VariableNames', varnames);

%%
% Now calculate the fractions and convert them to percentages
quantified.frac_flight = quantified.n_flight ./ quantified.n_obs .* 100; % flight of all obs
quantified.frac_soar = quantified.n_soaring ./ quantified.n_obs .* 100; % soaring of all obs
quantified.frac_oroglift = quantified.n_oroglift ./ quantified.n_obs .* 100; % oroglift soaring of all obs
quantified.frac_soar_flight = quantified.n_soaring ./ quantified.n_flight .* 100; % soaring of all flight obs
quantified.frac_oroglift_flight = quantified.n_oroglift ./ quantified.n_flight .* 100; % oroglift soaring of all flight obs
quantified.frac_oroglift_soaring = quantified.n_oroglift ./ quantified.n_soaring .* 100; % oroglift soaring of all soaring obs
quantified.frac_oroglift_all_altitudes = quantified.n_oroglift_all ./ quantified.n_obs .* 100; % altitude independent orographic lift soaring of all obs
quantified.frac_oroglift_flight_all_altitudes = quantified.n_oroglift_all ./ quantified.n_flight .* 100; % altitude independent oroglift soaring of all flight obs
quantified.frac_oroglift_soaring_all_altitudes = quantified.n_oroglift_all ./ quantified.n_soaring .* 100; % altitude independent oroglift soaring of all soaring obs

% The following are in the thesis and prefixed with r_
quantified.r_frac_flight = quantified.n_flight ./ quantified.n_obs .* 100; % percentage flight of all observations
quantified.r_frac_soaring_flight = quantified.n_soaring ./ quantified.n_flight .* 100; % percentage soaring flight of all flight observations
quantified.r_frac_oroglift_all_altitudes = quantified.n_oroglift_all ./ quantified.n_soaring .* 100; % percentage soaring flight over landscape causing threshold values of orographic lift
quantified.r_frac_oroglift_soaring_threshold = quantified.n_oroglift ./ quantified.n_soaring .* 100; % percentage soaring flight over landscape causing threshold values of orographic lift below max_meters altitude

%% Night summary
quantified_night.r_frac_flight = quantified_night.n_flight ./ quantified_night.n_obs .* 100; % percentage flight of all observations
quantified_night.r_frac_soaring_flight = quantified_night.n_soaring ./ quantified_night.n_flight .* 100; % percentage soaring flight of all flight observations
quantified_night.r_frac_oroglift_all_altitudes = quantified_night.n_oroglift_all ./ quantified_night.n_soaring .* 100; % percentage soaring flight over landscape causing threshold values of orographic lift
quantified_night.r_frac_oroglift_soaring_threshold = quantified_night.n_oroglift ./ quantified_night.n_soaring .* 100; % percentage soaring flight over landscape causing threshold values of orographic lift below max_meters altitude

%% Plot some of the results
x = 1:1:max_meters;
y = mean(quantified.frac_oroglift_soaring);
err = std(quantified.frac_oroglift_soaring);

errorbar(x, y, err);

%% Lesser Black-backed Gulls
% Flight
% lbbg = quantified(strcmp(quantified.project, 'LBBG_TEXEL') == 1, :);
% m = mean(lbbg.frac_flight)
% s = std(lbbg.frac_flight)
% cv = s / m
% minimum = min(lbbg.frac_flight)
% maximum = max(lbbg.frac_flight)

% Soaring
% m = mean(lbbg.r_frac_soaring_flight)
% s = std(lbbg.r_frac_soaring_flight)
% cv = s / m
% minimum = min(lbbg.r_frac_soaring_flight)
% maximum = max(lbbg.r_frac_soaring_flight)

% Soaring over landscape causing orographic lift
% m = mean(lbbg.r_frac_oroglift_all_altitudes)
% s = std(lbbg.r_frac_oroglift_all_altitudes)
% cv = s / m
% minimum = min(lbbg.r_frac_oroglift_all_altitudes)
% maximum = max(lbbg.r_frac_oroglift_all_altitudes)

% Soaring over landscape causing orographic lift at an altitude up to 25m
% l25 = lbbg.r_frac_oroglift_soaring_threshold(:, 25);
% m = mean(l25)
% s = std(l25)
% cv = s / m
% minimum = min(l25)
% maximum = max(l25)

% Soaring over landscape causing orographic lift at an altitude up to 10m
% l10 = lbbg.r_frac_oroglift_soaring_threshold(:, 10);
% m = mean(l10)
% s = std(l10)
% cv = s / m
% minimum = min(l10)
% maximum = max(l10)

%% Herring Gulls
% Flight
% hg = quantified(strcmp(quantified.project, 'HG_TEXEL') == 1, :);
% m = mean(hg.frac_flight)
% s = std(hg.frac_flight)
% cv = s / m
% minimum = min(hg.frac_flight)
% maximum = max(hg.frac_flight)

% Soaring
% m = mean(hg.r_frac_soaring_flight)
% s = std(hg.r_frac_soaring_flight)
% cv = s / m
% minimum = min(hg.r_frac_soaring_flight)
% maximum = max(hg.r_frac_soaring_flight)

% Soaring over landscape causing orographic lift
% m = mean(hg.r_frac_oroglift_all_altitudes)
% s = std(hg.r_frac_oroglift_all_altitudes)
% cv = s / m
% minimum = min(hg.r_frac_oroglift_all_altitudes)
% maximum = max(hg.r_frac_oroglift_all_altitudes)

% Soaring over landscape causing orographic lift at an altitude up to 25m
% l25 = hg.r_frac_oroglift_soaring_threshold(:, 25);
% m = mean(l25)
% s = std(l25)
% cv = s / m
% minimum = min(l25)
% maximum = max(l25)

% Soaring over landscape causing orographic lift at an altitude up to 10m
% l10 = hg.r_frac_oroglift_soaring_threshold(:, 10);
% m = mean(l10)
% s = std(l10)
% cv = s / m
% minimum = min(l10)
% maximum = max(l10)

%% Compare species
ranksum(lbbg.r_frac_flight, hg.r_frac_flight); % p = 0.9871 - No significant difference between groups
ranksum(lbbg.r_frac_soaring_flight, hg.r_frac_soaring_flight); % p = 0.717 - No significant difference between groups
ranksum(lbbg.r_frac_oroglift_all_altitudes, hg.r_frac_oroglift_all_altitudes); % p = 0.0033 - Significant difference between groups

ql25hg = hg.r_frac_oroglift_soaring_threshold(:, 25);
ql25lbbg = lbbg.r_frac_oroglift_soaring_threshold(:, 25);
ranksum(ql25hg, ql25lbbg); % p = 0.0055 - Significant differences between groups

ql10hg = hg.r_frac_oroglift_soaring_threshold(:, 10);
ql10lbbg = lbbg.r_frac_oroglift_soaring_threshold(:, 10);
ranksum(ql10hg, ql10lbbg); % p = 0.0099 - Significant differences between groups

%% Night summary statistics: Both species
qn = quantified_night;

% Flight by night
m = nanmean(qn.r_frac_flight)
s = nanstd(qn.r_frac_flight)
cv = s / m
minimum = nanmin(qn.r_frac_flight)
maximum = nanmax(qn.r_frac_flight)

% Soaring by night
m = nanmean(qn.r_frac_soaring_flight)
s = nanstd(qn.r_frac_soaring_flight)
cv = s / m
minimum = nanmin(qn.r_frac_soaring_flight)
maximum = nanmax(qn.r_frac_soaring_flight)

% Soaring over landscape causing orographic lift by night
m = nanmean(qn.r_frac_oroglift_all_altitudes)
s = nanstd(qn.r_frac_oroglift_all_altitudes)
cv = s / m
minimum = nanmin(qn.r_frac_oroglift_all_altitudes)
maximum = nanmax(qn.r_frac_oroglift_all_altitudes)

% Soaring over landscape causing orographic lift at an altitude up to 25m
% by night
l25 = qn.r_frac_oroglift_soaring_threshold(:, 25);
m = nanmean(l25)
s = nanstd(l25)
cv = s / m
minimum = nanmin(l25)
maximum = nanmax(l25)

% Soaring over landscape causing orographic lift at an altitude up to 10m
% by night
l10 = qn.r_frac_oroglift_soaring_threshold(:, 10);
m = nanmean(l10)
s = nanstd(l10)
cv = s / m
minimum = nanmin(l10)
maximum = nanmax(l10)

%% Night Summary Statistics: Herring Gull
hg_night = quantified_night(strcmp(quantified_night.project, 'HG_TEXEL') == 1, :);

% Flight by night
m = nanmean(hg_night.r_frac_flight)
s = nanstd(hg_night.r_frac_flight)
cv = s / m
minimum = nanmin(hg_night.r_frac_flight)
maximum = nanmax(hg_night.r_frac_flight)

% Soaring by night
m = nanmean(hg_night.r_frac_soaring_flight)
s = nanstd(hg_night.r_frac_soaring_flight)
cv = s / m
minimum = nanmin(hg_night.r_frac_soaring_flight)
maximum = nanmax(hg_night.r_frac_soaring_flight)

% Soaring over landscape causing orographic lift by night
m = nanmean(hg_night.r_frac_oroglift_all_altitudes)
s = nanstd(hg_night.r_frac_oroglift_all_altitudes)
cv = s / m
minimum = nanmin(hg_night.r_frac_oroglift_all_altitudes)
maximum = nanmax(hg_night.r_frac_oroglift_all_altitudes)

% Soaring over landscape causing orographic lift at an altitude up to 25m
% by night
l25 = hg_night.r_frac_oroglift_soaring_threshold(:, 25);
m = nanmean(l25)
s = nanstd(l25)
cv = s / m
minimum = nanmin(l25)
maximum = nanmax(l25)

% Soaring over landscape causing orographic lift at an altitude up to 10m
% by night
l10 = hg_night.r_frac_oroglift_soaring_threshold(:, 10);
m = nanmean(l10)
s = nanstd(l10)
cv = s / m
minimum = nanmin(l10)
maximum = nanmax(l10)

%% Night Summary Statistics: Lesser Black-backed Gull
lbbg_night = quantified_night(strcmp(quantified_night.project, 'LBBG_TEXEL') == 1, :);

% Flight by night
m = nanmean(lbbg_night.r_frac_flight)
s = nanstd(lbbg_night.r_frac_flight)
cv = s / m
minimum = nanmin(lbbg_night.r_frac_flight)
maximum = nanmax(lbbg_night.r_frac_flight)

% Soaring by night
m = nanmean(lbbg_night.r_frac_soaring_flight)
s = nanstd(lbbg_night.r_frac_soaring_flight)
cv = s / m
minimum = nanmin(lbbg_night.r_frac_soaring_flight)
maximum = nanmax(lbbg_night.r_frac_soaring_flight)

% Soaring over landscape causing orographic lift by night
m = nanmean(lbbg_night.r_frac_oroglift_all_altitudes)
s = nanstd(lbbg_night.r_frac_oroglift_all_altitudes)
cv = s / m
minimum = nanmin(lbbg_night.r_frac_oroglift_all_altitudes)
maximum = nanmax(lbbg_night.r_frac_oroglift_all_altitudes)

% Soaring over landscape causing orographic lift at an altitude up to 25m
% by night
l25 = lbbg_night.r_frac_oroglift_soaring_threshold(:, 25);
m = nanmean(l25)
s = nanstd(l25)
cv = s / m
minimum = nanmin(l25)
maximum = nanmax(l25)

% Soaring over landscape causing orographic lift at an altitude up to 10m
% by night
l10 = lbbg_night.r_frac_oroglift_soaring_threshold(:, 10);
m = nanmean(l10)
s = nanstd(l10)
cv = s / m
minimum = nanmin(l10)
maximum = nanmax(l10)
%% Test
ranksum(lbbg_night.r_frac_flight, hg_night.r_frac_flight); % p = 0.9871 - No significant difference between groups
ranksum(lbbg_night.r_frac_soaring_flight, hg_night.r_frac_soaring_flight); % p = 0.717 - No significant difference between groups
ranksum(lbbg_night.r_frac_oroglift_all_altitudes, hg_night.r_frac_oroglift_all_altitudes); % p = 0.0033 - Significant difference between groups

ql25hg = hg_night.r_frac_oroglift_soaring_threshold(:, 25);
ql25lbbg = lbbg_night.r_frac_oroglift_soaring_threshold(:, 25);
ranksum(ql25hg, ql25lbbg); % p = 0.0055 - Significant differences between groups

ql10hg = hg_night.r_frac_oroglift_soaring_threshold(:, 10);
ql10lbbg = lbbg_night.r_frac_oroglift_soaring_threshold(:, 10);
ranksum(ql10hg, ql10lbbg); % p = 0.0099 - Significant differences between groups

test(