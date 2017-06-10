clear;
close all;
clc;
%%
addpath(genpath('functions/'));

tracks = load('data/flight_tracks_2015_2016.mat');
tracks = tracks.flight_tracks;

% Make sure the table is sorted
tracks = sortrows(tracks,{'year', 'ID'},{'ascend','ascend'});

%% Add altitude above DEM level column to dataset
% Since we assume orographic lift soaring happens close to the surface,
% this is a relevant column to use to select the right datapoints.
tracks.altitude_adl = tracks.altitude - tracks.dem_alt_mean;

%% Summarize tracks
devices = unique(tracks.device_info_serial);
n_devices = size(devices, 1);
years = unique(tracks.year);
n_years = size(years, 1);

% Orographic lift parameters
max_meters = 30;
oroglift_min = 0.9;
oroglift_max = 4.2;

% Prepare columns
n_obs_year = [];
n_flight_obs_year = [];
n_soaring_obs_year = [];
n_orog_soaring_all_year = [];
n_orog_soaring_obs_year = zeros(1, max_meters);
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
        
        % Select soaring
        soaring_sel = flight_sel(flight_sel.class_id == 3, :);
        n_soaring_obs_year = [n_soaring_obs_year; size(soaring_sel, 1)];
        
        % Select orographic soaring without altitude
        % Criteria:
        % 1. Oroglift over 0.9m/s and below 4.2m/s
        orog_soaring_sel = soaring_sel(soaring_sel.oroglift_max >= oroglift_min & ...
                           soaring_sel.oroglift_max <= oroglift_max, :);
        n_orog_soaring_all_year = [n_orog_soaring_all_year; size(orog_soaring_sel, 1)];
        
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
        
    end
end

n_orog_soaring_obs_year = n_orog_soaring_obs_year(2:end, :);

varnames = {'device', 'year', 'project', 'n_obs', 'n_flight', ...
            'n_soaring', 'n_oroglift_all', 'n_oroglift'};
        
quantified = table(device, year, project, n_obs_year, n_flight_obs_year, ...
                   n_soaring_obs_year, n_orog_soaring_all_year, n_orog_soaring_obs_year, ...
                   'VariableNames', varnames);

% Now calculate the fractions and convert them to percentages
quantified.frac_flight = quantified.n_flight ./ quantified.n_obs .* 100; % flight of all obs
quantified.frac_soar = quantified.n_soaring ./ quantified.n_obs .* 100; % soaring of all obs
quantified.frac_soar_flight = quantified.n_soaring ./ quantified.n_flight .* 100; % soaring of all flight obs
quantified.frac_oroglift = quantified.n_oroglift ./ quantified.n_obs .* 100; % oroglift soaring of all obs
quantified.frac_oroglift_flight = quantified.n_oroglift ./ quantified.n_flight .* 100; % oroglift soaring of all flight obs
quantified.frac_oroglift_soaring = quantified.n_oroglift ./ quantified.n_soaring .* 100; % oroglift soaring of all soaring obs
quantified.frac_oroglift_all = quantified.n_oroglift_all ./ quantified.n_obs .* 100; % altitude independent orographic lift soaring of all obs
quantified.frac_oroglift_flight = quantified.n_oroglift_all ./ quantified.n_flight .* 100; % altitude independent oroglift soaring of all flight obs
quantified.frac_oroglift_soaring_all = quantified.n_oroglift_all ./ quantified.n_soaring .* 100; % altitude independent oroglift soaring of all soaring obs


%% Plot some of the results
x = 1:1:max_meters;
y = mean(quantified.frac_oroglift_soaring);
err = std(quantified.frac_oroglift_soaring);

errorbar(x, y, err);


