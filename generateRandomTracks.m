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