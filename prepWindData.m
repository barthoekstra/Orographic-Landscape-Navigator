function [data, stations] = prepWindData(storpath, extent, datestart, dateend)
%prepWeatherData Downloads, cleans and stores KNMI weather data
%   Inputs:
%   1. A path to store the used files at (e.g. data/weather/)
%   2. A shapefile containing the extent of the research area, which will
%      be used to select the right weather stations
%   3. Start date in yyyymmdd format
%   4. End date in yyyymmdd format
%
%   Outputs:
%   1. A table containing weather data from selected weather stations in 
%      the given timeframe.
%   2. A table containing information about the selected weather stations
%
%   Example usage:
%   [data, stations] = prepWeatherData('data/weather/', ...
%                           'data/ResearchAreaNH.shp', 20140501, 20140715);
%
%   Weather data comes from the Dutch Meteorological Institute (KNMI). The
%   script below downloads a .nc file from the following url containing all
%   prerequisite information about the weather stations. We need to use
%   this file to know the GPS locations of the weather stations and to
%   potentially use weather data from weather stations located at sea
%   (unfortunately not possible yet)
%
%   Alternative download if the script does not work:
%   1) Go to the following url to find the NetCDF file:
%      https://data.overheid.nl/data/dataset/a4c85871-6847-4b55-bc46-e14d74cf2c05
%      You should be on the KNMI network of observation stations dataset
%      page. If the URL does not work anymore, for some reason, search for
%      'KNMI network of observation stations' using the search function on
%      data.overheid.nl
%   2) Open the ATOM feed. 
%   3) Find the url that ends in .nc.zip and use it instead of the one 
%      below if it does not work
%
%   The following scripts automates these steps, but these steps may change
%   because the url may change. Tthe steps taken, however,  are not very
%   likely to change significantly.

    % ---------------------------------------------------------------------
    % Obtain NetCDF file containing information about all the automatic
    % weather stations at our disposal
    %
    % Note: although this file contains data about off-shore weather
    % stations, stations, the script is not able to download the data of
    % these stations subsequently yet. Maybe a work-around can be found,
    % so let's ignore these stations for now.
    % ---------------------------------------------------------------------

    if exist(storpath, 'dir') ~= 7
        mkdir(storpath);
    end

    if exist([storpath, 'obsStations.nc'], 'file') ~= 2
        obss = 'http://data.knmi.nl/inspire/download/waarneemstations/2/noversion/0000/00/00/STD___OPER_P___OBS_____L2.nc.zip';
        filename = 'obsStations.nc.zip';
        zipped = websave([storpath, filename], obss);
        unzipped = unzip(zipped, storpath);
        movefile(char(unzipped),[storpath, 'obsStations.nc']); % rename
        delete(zipped);
    end
    
    % ---------------------------------------------------------------------
    % Read NetCDF and select the weather stations that fall inside the
    % shapefile which demarcates our research area
    % ---------------------------------------------------------------------

    % Note: It turns out MATLAB does not support the newer implementation
    % of NetCDF yet, so we need to resort to using the HDF5 functions,
    % which functions identical.

    name = h5read([storpath,'obsStations.nc'], '/name/');  % station name
    wmo = h5read([storpath,'obsStations.nc'], '/WMO/');    % station ID
    type = h5read([storpath,'obsStations.nc'], '/type/');  % station type
    lat = h5read([storpath,'obsStations.nc'], '/lat/');    % station lat
    lon = h5read([storpath,'obsStations.nc'], '/lon/');    % station lon

    weatherStationsAll = table(wmo, name, type, lat, lon, 'RowNames', wmo);

    % Adjust wmo (weather station IDs) to fit use for later on by removing
    % 06-prefix
    weatherStationsAll(:,1) = cellfun(@(s) s(3:end), ...
        table2cell(weatherStationsAll(:,1)), 'uniformoutput', 0);
    weatherStationsAll.Properties.RowNames = weatherStationsAll.wmo;

    % Select only weather stations within research area polygon bounds
    polygon = shaperead(extent);

    in = inpolygon(lon, lat, polygon.X, polygon.Y);
    weatherStations = weatherStationsAll(in, :);

    % ---------------------------------------------------------------------
    % Download weather data from the KNMI using the selected weather 
    % stations as input
    % ---------------------------------------------------------------------
    
    % See information about usage of the following request method here:
    % http://www.knmi.nl/kennis-en-datacentrum/achtergrond/data-ophalen-vanuit-een-script
    url = 'http://projects.knmi.nl/klimatologie/uurgegevens/getdata_uur.cgi';

    % The NetCDF datasets contains station IDs which we need, but prefixes
    % them with '06', which we need to take off, because base_url used them
    % without this prefix. Also stations have to be separated with a colon
    % (:)
    stations = strjoin(weatherStations.wmo, ':');

    % Prepare other POST request components
    vars = 'WIND';                  % Variables of interest (see docs)
    datestartstr = [num2str(datestart), '01'];   % Start date with hour
    dateendstr   = [num2str(dateend), '24'];     % End date with hour

    % POST request example:
    % stns=235:280:260&vars=VICL:PRCP&start=2014050101&end=2014071524

    % Download file
    filename = [storpath, 'winddata.txt'];
    options = weboptions('RequestMethod', 'post');
    try
        weatherdata = websave(filename, url, 'stns', stations, 'vars', ...
                        vars, 'start', datestartstr, 'end', dateendstr, options);
    catch E
        disp(E.message);
    end
    
    % ---------------------------------------------------------------------
    % Data has been downloaded, now we can read it and remove the empty
    % data from the dataset
    % ---------------------------------------------------------------------
    
    % KNMI weather data always contains quite a few headerlines. In case
    % this script is only used to download wind data, the following is
    % correct. If more data is selected (by changing vars (above), the
    % following lines need to change
    headerlines = 15 + size(weatherStations,1);

    % Scan the contents of the file
    fileID = fopen(filename);
    weather = textscan(fileID, '%f%f%f%f%f%f%f', 'Delimiter', ',', 'EmptyValue',...
                       NaN, 'HeaderLines', headerlines);
    fclose(fileID);

    weather = cell2mat(weather);

    % Apparently this script can sometimes result in empty lines, as the
    % request to the server still generates a result for weather stations
    % with no data during the given timeframe. So we have to remove these
    % lines entirely. Most likely this will remove ALL the records of that
    % particular weather station, but maybe that's not always the case
    weatherStationsOld = unique(weather(:,1));
    weather = weather(~isnan(weather(:,4)),:);
    weatherStationsNew = unique(weather(:,1));

    % If ALL data from a weather station is unavailable, remove this
    % weather station from the overview of weather stations within the
    % polygon as well
    weatherStationDiff = setdiff(weatherStationsOld, weatherStationsNew);

    weatherStations.wmo = str2double(weatherStations.wmo);
    unusedWeatherStation = weatherStations.wmo == weatherStationDiff;
    weatherStations(unusedWeatherStation,:) = [];
    
    % Return values
    stations = weatherStations;
    data = table(weather(:,1), weather(:,2), weather(:,3), weather(:,4), weather(:,5), weather(:,6), weather(:,7), ...
             'VariableNames', {'stationID', 'date', 'hour', 'wdir', 'wspeed_hr', 'wspeed_10min', 'wspeed_peak'});
end

