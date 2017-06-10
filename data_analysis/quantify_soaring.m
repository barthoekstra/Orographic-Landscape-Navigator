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

%% Count consecutive number of observations within a day
% Essentially a run-length encoding and decoding problems, see function
% rude for documentation.
date = datenum(tracks.date_time, 'yyyy-mm-dd HH:MM:SS.0');
date = floor(date);
[date_consec, ~] = rude(date);
date_consec = rude(date_consec, date_consec);

% Now add it to the table
tracks.flight_obs_day = date_consec';

%% Calculate flight fraction
% Add flight fraction as a share of thumber of observations per day to the
% table
% tracks.flight_fraction = tracks.flight_obs_day ./ tracks.n_obs_day;
% 
% % Because the trackers do not log data at a constant rate and equally for
% % all individuals, we have to aggregate the values of all the days over all
% % the individuals in order to quantify soaring. First we will iterate over
% % the entire flight tracks dataset and summarize every day
% devices = unique(tracks.device_info_serial);
% years = unique(tracks.year);
% 
% device_col = [];
% year_col = [];
% month_col = [];
% day_col = [];
% flight_obs_day_col = [];
% soaring_obs_day_col = [];
% n_obs_day_col = [];
% 
% % Since we iterate over every day, the stepsize is equal to the amount of
% % flight observations in a single day (all other observations are excluded
% % anyways)
% [~, steps] = rude(tracks.flight_obs_day);
% 
% n = size(tracks, 1);
% k = 1;
% i = 1;
% 
% while i <= n & k <= size(steps, 2)
%     if steps(k) > 1
%         selection = tracks(i:(i + steps(k) - 1), :);
%     else
%         selection = tracks(i, :);
%         disp('Ding');
%     end
%     
%     device_col(k) = selection(1,:).device_info_serial;
%     year_col(k) = selection(1,:).year;
%     month_col(k) = selection(1,:).month;
%     day_col(k) = selection(1,:).day;
%     n_obs_day_col(k) = selection(1,:).n_obs_day;
%     flight_obs_day_col(k) = selection(1,:).flight_obs_day;
%     soaring_obs_day_col(k) = size(selection.class_id(selection.class_id == 3), 1);
%     
%     i = i + steps(k);
%     k = k + 1;
% end
% 
% % Prepare summary table
% varnames = {'device_info_serial', 'year', 'month', 'day', 'n_obs_day', ...
%             'flight_obs_day', 'soaring_obs_day'};
% 
% summary = table(device_col', year_col', month_col', day_col', ...
%                 n_obs_day_col', flight_obs_day_col', soaring_obs_day_col', ...
%                 'VariableNames', varnames);
% 
% % Add fractions
% summary.flight_fraction = summary.flight_obs_day ./ summary.n_obs_day;
% summary.soaring_fraction = summary.soaring_obs_day ./ summary.n_obs_day;
% summary.soaring_flight_fraction = summary.soaring_obs_day ./ summary.flight_obs_day;
% 
% % And summarise this table per individual per year
% statarray = grpstats(summary, {'device_info_serial', 'year'}, ... % group by
%                      {'min', 'max', 'mean', 'std'},...            % stats choice
%                      'DataVars', {'flight_fraction', 'soaring_fraction', 'soaring_flight_fraction'});

%% Summarize tracks
devices = unique(tracks.device_info_serial);
n_devices = size(devices, 1);
years = unique(tracks.year);
n_years = size(years, 1);

% Prepare columns
n_obs_year = [];
n_flight_obs_year = [];
n_soaring_obs_year = [];
n_orog_soaring_obs_year = [];
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
        
        % Select orographic soaring
        % Criteria:
        % 1. Oroglift over 0.9 m/s and below 4.2 m/s
        % 2. Altitude_adl <= 25
        orog_soaring_sel = soaring_sel(soaring_sel.oroglift_max >= 0.9 & ...
                                       soaring_sel.oroglift_max <= 4.2 & ...
                                       soaring_sel.altitude_adl, :);
        n_orog_soaring_obs_year = [n_orog_soaring_obs_year; size(orog_soaring_sel, 1)];
    end
end

varnames = {'device', 'year', 'n_obs', 'n_flight', ...
            'n_soaring', 'n_oroglift', 'project'};
        
quantified = table(device, year, n_obs_year, n_flight_obs_year, ...
                   n_soaring_obs_year, n_orog_soaring_obs_year, project, ...
                   'VariableNames', varnames);
                   
% Now calculate the fractions and convert them to percentages
quantified.frac_flight = quantified.n_flight ./ quantified.n_obs .* 100; % flight of all obs
quantified.frac_soar = quantified.n_soaring ./ quantified.n_obs .* 100; % soaring of all obs
quantified.frac_soar_flight = quantified.n_soaring ./ quantified.n_flight .* 100; % soaring of all flight obs
quantified.frac_oroglift = quantified.n_oroglift ./ quantified.n_obs .* 100; % oroglift soaring of all obs
quantified.frac_oroglift_soaring = quantified.n_oroglift ./ quantified.n_soaring .* 100; % oroglift soaring of all soaring obs

%% Alternative implementation of what is above
% Apparently the grpstats function can do what we have written above and it
% does so much faster.
% init_group = grpstats(tracks, ...                                                       % data source
%                       {'device_info_serial', 'year', 'month', 'day', 'class_id', 'n_obs_day'}, ...   % group by
%                       {'numel'}, ...                                                    % statistics
%                       'DataVars', {'class_id'});                           % 
% 
% % Rename column and remove unused columns
% init_group.n_obs_class = init_group.GroupCount;
% init_group.GroupCount = [];
% init_group.numel_class_id = [];
% 
% % Calculate fractions
% init_group.frac_class = init_group.n_obs_class ./ init_group.n_obs_day;
% 
% % Now aggregate this to a higher level: per bird per year per class type
% agg_group_classes = grpstats(init_group, ...
%                      {'device_info_serial', 'year', 'class_id'}, ...
%                      {'mean', 'min', 'max'}, ...
%                      'DataVars', {'frac_class'});
                  
%%
% Remove repetitions, again using rude, explained:
% The table now contains 3 times the flight_fraction for a day that has a
% flight_obs_day of 3. In order to analyse the fraction, we should thus
% remove this fraction 11 times and work with the remainder. So:
% [3 3 3 4 4 4 4] should become [3 4]
%[~, flight_fraction] = rude(tracks.flight_fraction);



% Because 