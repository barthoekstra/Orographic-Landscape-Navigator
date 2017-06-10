% In order to determine the threshold levels for orographic lift soaring,
% we want to create flight tracks with entirely random lat and lon
% coordinates. Following this approach, we can see what the difference is
% between the background orographic lift levels and the ones the birds
% choose for flight. Randomizing the lat and lon values of the tracks would
% represent behaviour as if the birds do not choose their flight paths with
% a certain strategy in mind. To keep all things equal (wind directions and
% speeds), we only randomize the lat and lon values.

load('proj_settings_2016.mat');
load('proj_tracks_2016.mat');

%%
% We want to generate lat and lon values within the same bounding box as
% the research area.
area = shapeinfo('data/ResearchAreaNH.shp');
bb = area.BoundingBox;
lon_min = bb(1,1);
lon_max = bb(2,1);
lat_min = bb(1,2);
lat_max = bb(2,2);

% Now generate uniformly distributed random numbers with these values for
% every datapoint in the tracks matrix and replace the values of the
% corresponding columns.
n = size(tracks, 1);
random_lon = (lon_max - lon_min) .* rand(n, 1) + lon_min;
random_lat = (lat_max - lat_min) .* rand(n, 1) + lat_min;

tracks.longitude = random_lon;
tracks.latitude = random_lat;

save(['proj_tracks_random_', daterange, '.mat'], 'tracks');

% This will generate a lot of points over sea - obviously - but for now
% this is the quickest and easiest approach, even if it takes longer for
% the computer to process the results. Points over the sea will be excluded
% later on anyways.

profile on;

load('proj_settings_2016.mat');
load('proj_tracks_random_2016.mat');
load('proj_dems_2016.mat');
load('proj_wind_2016.mat');


% Prepare tracks to store orographic lift data
i = size(tracks, 1);
tracks.oroglift_max = zeros(i,1);
tracks.oroglift_mean = zeros(i,1);
tracks.oroglift_min = zeros(i,1);
tracks.dem_alt_max = zeros(i,1);
tracks.dem_alt_mean = zeros(i,1);
tracks.dem_alt_min = zeros(i,1);
tracks.wspeed = zeros(i,1);
tracks.wdir = zeros(i,1);
tracks.wstation = zeros(i,1);
tracks.asp_mode = zeros(i,1);
tracks.asp_std = zeros(i,1);

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

save(['proj_tracks_oroglift_random_', daterange, '.mat'], 'tracks');

profile viewer;

%% Merge tracks and select random points
% Now that we have generated random points with corresponding orographic
% lift values, we will select a random subset of both years (2015 and 2016)
% which are equal in length to the true datasets of both years. Count the
% number of flight observations of both years and put them in the variables
% below:
flight_tracks = load('flight_tracks_2015_2016.mat');
flight_tracks = flight_tracks.flight_tracks;

nr2015 = sum(flight_tracks.year == 2015);
nr2016 = sum(flight_tracks.year == 2016);

% Now we select nr2015 random datapoints from the 2015 dataset and nr2016
% from the 2016 dataset
random2015 = load('proj_tracks_oroglift_random_2015.mat');
random2015 = random2015.tracks;
random2016 = load('proj_tracks_oroglift_random_2016.mat');
random2016 = random2016.tracks;

size2015 = size(random2015, 1);
size2016 = size(random2016, 1);

random2015 = random2015(randperm(size2015, nr2015), :);
random2016 = random2016(randperm(size2016, nr2016), :);

random_tracks = vertcat(random2015, random2016);
save('random_tracks_2015_2016', 'random_tracks');
writetable(random_tracks, 'random_tracks_2015_2016.csv');
