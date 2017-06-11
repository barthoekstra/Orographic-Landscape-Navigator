clear;
close all;
clc;

addpath(genpath('functions/'));

tracks = load('data/flight_tracks_2015_2016.mat');
tracks = tracks.flight_tracks;
random_tracks = load('data/random_tracks_2015_2016.mat');
random = random_tracks.random_tracks;
clear random_tracks;

%% Add altitude above DEM level column to dataset
% Since we assume orographic lift soaring happens close to the surface,
% this is a relevant column to use to select the right datapoints.
tracks.altitude_adl = tracks.altitude - tracks.dem_alt_mean;
random.altitude_adl = random.altitude - random.dem_alt_mean;

%% Remove datapoints at sea
% We are not interested in these datapoints, so we will remove them
% beforehand, so they don't affect the results. We can simply remove them
% by removing the points where the DEM altitude is exactly 0, without any
% decimals. This is so unlikely in reality that we can assume this only
% happens at sea, where we have set NaN values to 0.
tracks = tracks(tracks.dem_alt_max ~= 0, :);
random = random(random.dem_alt_max ~= 0, :);

%% Make flapping and soaring selections
% Classified soaring tracks have class 3 and classified flapping and
% extreme flapping tracks have class 1 and 2 respectively. Since manouvre,
% the last class type in this dataset is intermediate between soaring and
% flapping, we leave this out.
soaring = tracks(tracks.class_id == 3, :);
flapping = tracks(tracks.class_id == 1 | tracks.class_id == 2, :);

%% Determine consecutive track lengths
% In order to find high-resolution and highly accurate soaring tracks, we
% determine the length of consecutively increasing ID integers. So [1, 2,
% 3, 4, 6] result in [4, 1]. The following code calculates these values,
% which we can subsequently use to select long tracks of consecutive
% numbers.
IDs = soaring.ID;
consecutive_ints = diff(find(diff([nan; IDs(:); nan]) ~= 1));
consecutive_expanded = rude(consecutive_ints, consecutive_ints);
soaring.consecutive = consecutive_expanded';
soaring = [soaring(:,1), soaring(:,40), soaring(:,2:39)];

%% Select soaring datapoints
% In order to determine the characteristics of orographic lift soaring, we
% need to make a selection of high-resolution tracks of soaring flight,
% because those offer the most accurate data. We thus: 
% 1. Select soaring tracks with consecutively more than 5 datapoints which
%    are classified as soaring by the classifier.
% 2. Select only datapoints where the altitude above the mean DEM altitude
%    in meters is higher than adl_min and lower than adl_max, because we
%    assume orographic lift soaring happens fairly close to the surface
consecutive_min = 5;
adl_min = 0; % adl = above dem level
adl_max = 25;

hr_soaring = soaring(soaring.altitude_adl > adl_min & ...
                     soaring.altitude_adl <= adl_max & ...
                     soaring.consecutive >= consecutive_min, :);

%% Visualise orographic lift
% Now that we have loaded all data, including the randomized tracks, we can
% start comparing soaring conditions with those experienced by birds with
% no strategy with regards to flight at all (the random tracks). Let's
% first look at the difference between randomized flight, with a sort of
% background distribution of orographic lift, and strategic soaring flight
% by gulls

% First set some general visualisation settings
width = 6; % inch
height = 4; % inch
axislinewidth = 1; % 0.75 for paper
fontsize = 14; % 10 for paper
linewidth = 2; % 1.5 for paper
markersize = 12; % 8 for paper

% Now create the plot for both species combined
pd_soaring = fitdist(hr_soaring.oroglift_max, 'Kernel');
pd_random = fitdist(random.oroglift_max, 'Kernel');
x = -2:.1:10;
y_both_species = pdf(pd_soaring, x); % Both species of gulls combined
y_random = pdf(pd_random, x);

[x_hr_soaring_all, ~] = intersections(x, y_random, x, y_both_species, 1)

figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width * 100, height * 100]);
set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
hold on;
plot(x, y_random, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
plot(x, y_both_species, 'DisplayName', 'Actual tracks', 'LineWidth', 2);
title({'\bf\fontsize{14} Herring & Lesser Black-backed Gulls', ...
       '\rm\fontsize{12} Orographic lift soaring rates: 0.95 - 4.20 [m/s]'}); % add these values manually
xlabel('Orographic lift [m/s]', 'FontSize', fontsize); 
ylabel('Probability density', 'FontSize', fontsize);
legend('show');
hold off;

%% Let's now look at the differences between species
% Basically an identical approach to the previous section, but now we look
% at the species seperately.
hr_soaring_hg = hr_soaring(strcmp(hr_soaring.project, 'HG_TEXEL') == 1, :);
hr_soaring_lbbg = hr_soaring(strcmp(hr_soaring.project, 'LBBG_TEXEL') == 1, :);

pd_soaring_hg = fitdist(hr_soaring_hg.oroglift_max, 'Kernel');
pd_soaring_lbbg = fitdist(hr_soaring_lbbg.oroglift_max, 'Kernel');

y_hg = pdf(pd_soaring_hg, x);
y_lbbg = pdf(pd_soaring_lbbg, x);

[x_hr_soaring_hg, ~] = intersections(x, y_random, x, y_hg, 1)
[x_hr_soaring_lbbg, ~] = intersections(x, y_random, x, y_lbbg, 1)

figure(2);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width * 200, height * 100]);
% Herring Gull
subplot(1,2,1);
    set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
    hold on;
    plot(x, y_random, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
    plot(x, y_hg, 'DisplayName', 'Actual tracks', 'LineWidth', 2);
    title({'\bf\fontsize{14} Herring Gull', ...
           '\rm\fontsize{12} Orographic lift soaring rates: 0.93 - 4.18 [m/s]'}); % add these values manually
    xlabel('Orographic lift [m/s]', 'FontSize', fontsize); 
    ylabel('Probability density', 'FontSize', fontsize);
    legend('show');
    hold off;
    
% Lesser Black-backed Gull
subplot(1,2,2);
    set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
    hold on;
    plot(x, y_random, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
    plot(x, y_lbbg, 'DisplayName', 'Actual tracks', 'LineWidth', 2);
    title({'\bf\fontsize{14} Lesser Black-backed Gull', ...
           '\rm\fontsize{12} Orographic lift soaring rates: 0.99 - 4.33 [m/s]'}); % add these values manually
    xlabel('Orographic lift [m/s]', 'FontSize', fontsize); 
    ylabel('Probability density', 'FontSize', fontsize);
    legend('show');
    hold off;
    
% And finally plot them together in 1 figure
figure(3);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width * 100, height * 100]);
set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
hold on;
plot(x, y_random, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
plot(x, y_hg, 'DisplayName', 'Herring Gull', 'LineWidth', 2);
plot(x, y_lbbg, 'DisplayName', 'Lesser Black-backed Gull', 'LineWidth', 2);
title({'\bf\fontsize{14} Species comparison'});
xlabel('Orographic lift [m/s]', 'FontSize', fontsize); 
ylabel('Probability density', 'FontSize', fontsize);
legend('show');
hold off;

%% Now do the same, but with ALL soaring instead of only consecutive soaring
% We have now seen that there is only a slight difference in soaring
% preferences between Herring and Lesser Black-backed Gulls when we only
% analyse the consecutive stretches of more than 5 measurements of
% orographic lift soaring. It could be that differences are more pronounced
% when we look at the entire dataset, so let's do that. We simply
% copy-paste what is above.

% Select a subset of the soaring dataset, but now only select for altitude
soaring = soaring(soaring.altitude_adl > adl_min & ...
                  soaring.altitude_adl <= adl_max, :);

% Create the plot for both species combined
pd_soaring = fitdist(soaring.oroglift_max, 'Kernel');
x = -2:.1:10;
y_both_species = pdf(pd_soaring, x); % Both species of gulls combined

[x_soaring_all, ~] = intersections(x, y_random, x, y_both_species, 1)

figure(4);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width * 100, height * 100]);
set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
hold on;
plot(x, y_random, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
plot(x, y_both_species, 'DisplayName', 'Actual tracks', 'LineWidth', 2);
title({'\bf\fontsize{14} Herring & Lesser Black-backed Gulls', ...
       '\rm\fontsize{12} Orographic lift soaring rates: 0.82 - 4.27 [m/s]'}); % add these values manually
xlabel('Orographic lift [m/s]', 'FontSize', fontsize); 
ylabel('Probability density', 'FontSize', fontsize);
legend('show');
hold off;

%% And the differences between species
% Basically an identical approach to the previous section, but now we look
% at the species seperately.
soaring_all_hg = soaring(strcmp(soaring.project, 'HG_TEXEL') == 1, :);
soaring_all_lbbg = soaring(strcmp(soaring.project, 'LBBG_TEXEL') == 1, :);

pd_soaring_all_hg = fitdist(soaring_all_hg.oroglift_max, 'Kernel');
pd_soaring_all_lbbg = fitdist(soaring_all_lbbg.oroglift_max, 'Kernel');

y_all_hg = pdf(pd_soaring_all_hg, x);
y_all_lbbg = pdf(pd_soaring_all_lbbg, x);

[x_soaring_all_hg, ~] = intersections(x, y_random, x, y_all_hg, 1)
[x_soaring_all_lbbg, ~] = intersections(x, y_random, x, y_all_lbbg, 1)

figure(5);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width * 200, height * 100]);
% Herring Gull
subplot(1,2,1);
    set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
    hold on;
    plot(x, y_random, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
    plot(x, y_all_hg, 'DisplayName', 'Actual tracks', 'LineWidth', 2);
    title({'\bf\fontsize{14} Herring Gull', ...
           '\rm\fontsize{12} Orographic lift soaring rates: 0.80 - 4.23 [m/s]'}); % add these values manually
    xlabel('Orographic lift [m/s]', 'FontSize', fontsize); 
    ylabel('Probability density', 'FontSize', fontsize);
    legend('show');
    hold off;
    
% Lesser Black-backed Gull
subplot(1,2,2);
    set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
    hold on;
    plot(x, y_random, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
    plot(x, y_all_lbbg, 'DisplayName', 'Actual tracks', 'LineWidth', 2);
    title({'\bf\fontsize{14} Lesser Black-backed Gull', ...
           '\rm\fontsize{12} Orographic lift soaring rates: 0.85 - 6.50 [m/s]'}); % add these values manually
    xlabel('Orographic lift [m/s]', 'FontSize', fontsize); 
    ylabel('Probability density', 'FontSize', fontsize);
    legend('show');
    hold off;
    
% And finally plot them together in 1 figure
figure(6);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width * 100, height * 100]);
set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
hold on;
plot(x, y_random, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
plot(x, y_all_hg, 'DisplayName', 'Herring Gull', 'LineWidth', 2);
plot(x, y_all_lbbg, 'DisplayName', 'Lesser Black-backed Gull', 'LineWidth', 2);
title({'\bf\fontsize{14} Species comparison'});
xlabel('Orographic lift [m/s]', 'FontSize', fontsize); 
ylabel('Probability density', 'FontSize', fontsize);
legend('show');
hold off;

%% Determine some summary statistics for the thresholds
% Since true thresholds have already been determined, we only calculate
% median, mean, std and cv
mean_all = nanmean(soaring.oroglift_max)
median_all = nanmedian(soaring.oroglift_max)
std_all = nanstd(soaring.oroglift_max)
cv_all = std_all / mean_all

mean_hg_all = nanmean(soaring_all_hg.oroglift_max)
median_hg_all = nanmedian(soaring_all_hg.oroglift_max)
std_hg_all = nanstd(soaring_all_hg.oroglift_max)
cv_hg_all = std_hg_all / mean_hg_all

mean_lbbg_all = nanmean(soaring_all_lbbg.oroglift_max)
median_lbbg_all = nanmedian(soaring_all_lbbg.oroglift_max)
std_lbbg_all = nanstd(soaring_all_lbbg.oroglift_max)
cv_lbbg_all = std_lbbg_all / mean_lbbg_all

%% Comparison with randomized tracks
% Visually a difference is pretty clear. Let's test that statistially as
% well, just to be sure:
[h, p, ks2stat] = kstest2(random.oroglift_max, hr_soaring.oroglift_max) % h = 1; p = 0; ks2stat = 0.4575
[h, p, ks2stat] = kstest2(random.oroglift_max, soaring.oroglift_max) % h = 1; p = 0; ks2stat = 0.3720

%% Compare species
% Visually, species appear to respond differently to orographic lift, but
% are these effects significant? Again we use Wilcoxon-Mann-Whitney to
% compare the groups
ranksum(soaring_all_hg.oroglift_max, soaring_all_lbbg.oroglift_max) % p = 5.4527e-13 - Significant difference between groups
[h, p, ks2stat] = kstest2(soaring_all_hg.oroglift_max, soaring_all_lbbg.oroglift_max) % h = 1; p = 4.5028e-48; ks2stat = 0.0775
ranksum(hr_soaring_hg.oroglift_max, hr_soaring_lbbg.oroglift_max) % p = 0.0012 - Significant difference between groups
[h, p, ks2stat] = kstest2(hr_soaring_hg.oroglift_max, hr_soaring_lbbg.oroglift_max) % h = 1; p = 3.0698e-05; ks2stat = 0.0362

% This corroborates the visuals: both groups are significantly different in
% both the hr_soaring and overall soaring groups.

%% Wind speed and wind direction relative to landscape aspect
% We have looked at orographic lift rates, but we can do something similar
% to determine the thresholds for orographic lift soaring in terms of wind
% speed and wind direction relative to the landscape aspect. Assessing wind
% direction on its own is not very relevant (although could be done),
% because it does not say anything about orographic lift

% Create the plot for both species combined
pd_soaring_windspeed = fitdist(soaring.wspeed, 'Kernel');
pd_random_windspeed = fitdist(random.wspeed, 'Kernel');
x = 0:0.5:15;
y_both_species_windspeed = pdf(pd_soaring_windspeed, x); % Both species of gulls combined
y_random_windspeed = pdf(pd_random_windspeed, x);

[x_soaring_windspeed, ~] = intersections(x, y_random_windspeed, x, y_both_species_windspeed, 1)

figure(7);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width * 100, height * 100]);
set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
hold on;
plot(x, y_random_windspeed, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
plot(x, y_both_species_windspeed, 'DisplayName', 'Actual tracks', 'LineWidth', 2);
title({'\bf\fontsize{14} Herring & Lesser Black-backed Gulls', ...
       '\rm\fontsize{12} Wind speeds: 4 - 11 [m/s]'}); % add these values manually
xlabel('Wind speed [m/s]', 'FontSize', fontsize); 
ylabel('Probability density', 'FontSize', fontsize);
legend('show');
hold off;

% Species seperately
pd_soaring_hg_windspeed = fitdist(soaring_all_hg.wspeed, 'Kernel');
pd_soaring_lbbg_windspeed = fitdist(soaring_all_lbbg.wspeed, 'Kernel');

y_hg_windspeed = pdf(pd_soaring_hg_windspeed, x);
y_lbbg_windspeed = pdf(pd_soaring_lbbg_windspeed, x);

[x_soaring_hg_windspeed, ~] = intersections(x, y_random_windspeed, x, y_hg_windspeed, 1)
[x_soaring_lbbg_windspeed, ~] = intersections(x, y_random_windspeed, x, y_lbbg_windspeed, 1)

figure(8);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width * 200, height * 100]);
% Herring Gull
subplot(1,2,1);
    set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
    hold on;
    plot(x, y_random_windspeed, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
    plot(x, y_hg_windspeed, 'DisplayName', 'Actual tracks', 'LineWidth', 2);
    title({'\bf\fontsize{14} Herring Gull', ...
           '\rm\fontsize{12} Wind speeds: 4 - 11 [m/s]'}); % add these values manually
    xlabel('Wind speed [m/s]', 'FontSize', fontsize); 
    ylabel('Probability density', 'FontSize', fontsize);
    legend('show');
    hold off;
    
% Lesser Black-backed Gull
subplot(1,2,2);
    set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
    hold on;
    plot(x, y_random_windspeed, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
    plot(x, y_lbbg_windspeed, 'DisplayName', 'Actual tracks', 'LineWidth', 2);
    title({'\bf\fontsize{14} Lesser Black-backed Gull', ...
           '\rm\fontsize{12} Wind speeds: 5 - 10 [m/s]'}); % add these values manually
    xlabel('Wind speed [m/s]', 'FontSize', fontsize); 
    ylabel('Probability density', 'FontSize', fontsize);
    legend('show');
    hold off;
    
% And finally plot them together in 1 figure
figure(9);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width * 100, height * 100]);
set(gca, 'FontSize', fontsize, 'LineWidth', axislinewidth);
hold on;
plot(x, y_random_windspeed, 'DisplayName', 'Randomized tracks', 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
plot(x, y_hg_windspeed, 'DisplayName', 'Herring Gull', 'LineWidth', 2);
plot(x, y_lbbg_windspeed, 'DisplayName', 'Lesser Black-backed Gull', 'LineWidth', 2);
title({'\bf\fontsize{14} Species comparison'});
xlabel('Wind speed [m/s]', 'FontSize', fontsize); 
ylabel('Probability density', 'FontSize', fontsize);
legend('show');
hold off;
