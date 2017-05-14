%clear;
%close all;
%clc;

% Weather data comes from the Dutch Meteorological Institute. The following
% script downloads a .nc file from the following url containing all
% prerequisite information about the weather stations. We need to use this
% to later on access data from weather stations at sea, by default
% acquiring non-land-based data is not documented anywhere (at least I 
% could not find it). 
% 
% 1) Go to the following url to find the NetCDF file:
%    https://data.overheid.nl/data/dataset/a4c85871-6847-4b55-bc46-e14d74cf2c0
%    You should be on the ?KNMI network of observation stations? dataset
%    page
% 2) Open the ATOM feed.
% 3) Find the url that ends in .nc.zip and use it instead of the one below
%    if it does not work
%
% The following scripts automates these steps, but these steps may change
% because the url may change. However, the steps taken are not very likely
% to change significantly.

profile on;

storpath = 'data/weather/';

% ---------------------------------------------------------------------
% Obtain NetCDF file containing information about all the automatic 
% weather stations at our disposal
% ---------------------------------------------------------------------

if exist(storpath,'dir') ~= 7
    mkdir(storpath);
end

if exist([storpath,'obsStations.nc'],'file') ~= 2
    obss = 'http://data.knmi.nl/inspire/download/waarneemstations/2/noversion/0000/00/00/STD___OPER_P___OBS_____L2.nc.zip';
    filename = 'obsStations.nc.zip';
    zipped = websave([storpath,filename], obss);
    unzipped = unzip(zipped,storpath);
    movefile(char(unzipped),[storpath,'obsStations.nc']); % rename
    delete(zipped);
end

% ---------------------------------------------------------------------
% Read NetCDF and select the weather stations that fall inside the
% shapefile which demarcates our research area
% ---------------------------------------------------------------------

% Note: It turns out MATLAB does not support the newer implementation of
% NetCDF yet, so we need to resort to using the HDF5 functions, which
% functions identical.

name = h5read([storpath,'obsStations.nc'], '/name/');
wmo = h5read([storpath,'obsStations.nc'], '/WMO/');
type = h5read([storpath,'obsStations.nc'], '/type/');
lat = h5read([storpath,'obsStations.nc'], '/lat/');
lon = h5read([storpath,'obsStations.nc'], '/lon/');

weatherStationsAll = table(wmo, name, type, lat, lon, 'RowNames',wmo);

% Adjust wmo (weather station IDs) to fit use later on by removing
% 06-prefix
weatherStationsAll(:,1) = cellfun(@(s) s(3:end), table2cell(weatherStationsAll(:,1)), 'uniformoutput', 0);
weatherStationsAll.Properties.RowNames = weatherStationsAll.wmo;

% Select only weather stations within polygon bounds
polygon = shaperead('data/ResearchAreaNH.shp');

in = inpolygon(lon, lat, polygon.X, polygon.Y);
weatherStations = weatherStationsAll(in,:);

% ---------------------------------------------------------------------
% Download weather data from the KNMI using the selected weather 
% stations as input
% ---------------------------------------------------------------------

% The following url requires the actual request to be sent through POST:
%%

% See information about usage of the following request method here:
% http://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
url = 'http://projects.knmi.nl/klimatologie/uurgegevens/getdata_uur.cgi';

% The NetCDF datasets contains station IDs which we need, but prefixes them
% with '06', which we need to take off, because base_url used them without
% this prefix. Also stations have to be separated with a colon (:)
stations = strjoin(weatherStations.wmo,':');
stations = replace(stations,'06','');

% Prepare other POST request components
vars = 'WIND';

dateformat = 'yyyymmdd';
datestartstr = '2014050101';
dateendstr   = '2014071524';

datestart = datevec(datestartstr, dateformat);
dateend   = datevec(dateendstr, dateformat);

% POST request example:
% stns=235:280:260&vars=VICL:PRCP&byear=1970&bmonth=1&bday=1&eyear=2009&emonth=8&eday=18
                   
% Download file
filename = [storpath,'winddata.txt'];
options = weboptions('RequestMethod', 'post');
try
    weatherdata = websave(filename, url, 'stns', stations, 'vars', vars, 'start', ...
                      datestartstr, 'end', dateendstr, options);
catch E
    disp(E);
end
%%
% Data downloaded, let's read it from the file

headerlines = 15 + size(weatherStations,1); % Change if more variables are used

fileID = fopen(filename);
weather = textscan(fileID, '%f%f%f%f%f%f%f', 'Delimiter', ',', 'EmptyValue',...
                   NaN, 'HeaderLines', headerlines);
fclose(fileID);

weather = cell2mat(weather);

% Apparently this script can sometimes result in empty lines, as the
% request to the server still generates a result for weather stations with
% no data during the given timeframe. So we have to remove these lines
% entirely. Most likely this will remove ALL the records of that particular
% weather station.
weatherStationsOld = unique(weather(:,1));
weather = weather(~isnan(weather(:,4)),:);
weatherStationsNew = unique(weather(:,1));

% If ALL data from a weather station is unavailable, remove this weather
% station from the overview of weather stations within the polygon as well
weatherStationDiff = setdiff(weatherStationsOld, weatherStationsNew);

weatherStations.wmo = str2double(weatherStations.wmo);
unusedWeatherStation = weatherStations.wmo == weatherStationDiff;
weatherStations(unusedWeatherStation,:) = [];


profile viewer




