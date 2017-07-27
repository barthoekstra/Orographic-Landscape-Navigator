function [used_devices] = prepBirdAccelerometerData(db_user, db_password, device_numbers, starttime, stoptime, minobs, shapefile, storpath)
%prepBirdAccelerometerData Queries accelerometer data
%   Inputs:
%   1. Database username
%   2. Database password
%   3. A vector containing tracker device IDs
%   4. A starttime string in the format 'yyyy-mm-dd hh:mm:ss'
%   5. An endtime string in the format 'yyyy-mm-dd hh:mm:ss'
%   6. A minimum number of observations per day the query will select
%   7. Path to a shapefile delineating the research area bounds

%   Example usage:
%   tracks = prepBirdAccelerometerData2(username, password, devices, '2014-05-01 00:00:00', '2014-07-15 00:00:00', 3600, 'data/ResearchAreaNH.shp', 'data/unclassified');
%
%   This script queries the database for accelerometer data for the
%   classifier. In order to make sure a sufficient density of data points
%   per day are collected, it uses a minimum number of observations per
%   day. Only days with at least this number of observations (with
%   corresponding accelerometer measurements) are returned by the database.
%   Consequently, if you want to have ALL the data, set it to 0.

    profile on;
    
    % Batching
    % Although not always necessary, the code is written to allow for
    % batching the download. I have not checked to what extent I can push
    % the limits of this code until I encounter memory issues, but the
    % following variables can be changed to do so. More info on the url
    % below:
    % https://nl.mathworks.com/help/database/ug/preference-settings-for-large-data-import.html
    batchsize = 25000; % change also in line below:
    setdbprefs('FetchBatchSize', '25000'); % HAS to be written in quotes for some reason
    setdbprefs('FetchInBatches', 'yes');
    
    % Determine query boundaries based on the research area shapefile
    shape = shapeinfo(shapefile);
    % Boundingbox contains a 2x2 matrix with minimum and maximum lat and
    % lon values, as follows:
    %         Col 1     Col 2
    % Row 1   min(lon)  min(lat)
    % Row 2   max(lon)  max(lat)
    boundingbox = shape.BoundingBox;
    
    % Setup connection with database
    % Beware: In order to connect, you may have to be on the UvA network or
    % using a VPN.
    javaaddpath('drivers/postgresql-42.0.0.jre7.jar'); % required Driver
    conn = database('eecology', db_user, db_password, ...
        'org.postgresql.Driver', 'jdbc:postgresql://db.e-ecology.sara.nl:5432/eecology?ssl=true&sslfactory=org.postgresql.ssl.NonValidatingFactory&');
    setdbprefs('DataReturnFormat', 'table');
    
    % Storpath needs to exist to store the files
    if exist(storpath, 'dir') ~= 7
        mkdir(storpath);
    end

    % Generate neatly chuncked dates
    dates = generateDateStrings(starttime, stoptime);
    
    % For later use we store the vector of used devices
    devices = [];
    
    for i = 1:numel(device_numbers) % Iterate through the devices
        for j = 1:numel(dates)-1    % Iterate through date steps
            
            fprintf('Checking tracker %d for data on daterange %s till %s \n', device_numbers(i), dates(j), dates(j+1))
            % Turn datevalues into SQL-ready strings:
            start = strcat('''', dates(j), '''');
            stop = strcat('''', dates(j+1), '''');
            
            sql = strcat('WITH high_res_days AS (', ...
                    ' SELECT s.device_info_serial, s.date_time, ', ...
                        ' s.longitude, s.latitude, s.altitude, ', ...
                        ' s.speed_2d, ', ...
                        ' COUNT(*) OVER (PARTITION BY DATE(s.date_time)) AS n_obs_day ', ...
                    ' FROM gps.ee_tracking_speed_limited s', ...
                    ' WHERE ( s.date_time >= ', start,  ')', ...
                        ' AND ( s.date_time <= ', stop,    ')' , ...
                        ' AND ( s.device_info_serial = ', num2str(device_numbers(i)), ')', ...
                        ' AND ( s.longitude >= ', num2str(boundingbox(1,1)), ')', ...
                        ' AND ( s.longitude <= ', num2str(boundingbox(2,1)), ')', ...
                        ' AND ( s.latitude >= ', num2str(boundingbox(1,2)), ')', ...
                        ' AND ( s.latitude <= ', num2str(boundingbox(2,2)), ')', ...
                        ' AND s.userflag = 0', ...
                ' )', ...
                ' SELECT s.device_info_serial, s.date_time, ', ...
                    ' date_part(''year''::text, t.date_time) AS year, ', ...
                    ' date_part(''month''::text, t.date_time) AS month, ', ...
                    ' date_part(''day''::text, t.date_time) AS day, ', ...
                    ' date_part(''hour''::text, t.date_time) AS hour, ', ...
                    ' date_part(''minute''::text, t.date_time) AS minute, ', ...
                    ' date_part(''second''::text, t.date_time) AS second, ', ...
                    ' s.speed_2d AS speed, ', ...
                    ' s.longitude, s.latitude, s.altitude, ', ...
                    ' t.speed as tspeed, ', ...
                    ' a.index, (a.x_acceleration-d.x_o)/d.x_s AS x_cal, ', ...
                    ' (a.y_acceleration-d.y_o)/d.y_s AS y_cal, ', ...
                    ' (a.z_acceleration-d.z_o)/d.z_s AS z_cal ', ...
                ' FROM high_res_days s' , ...
                ' LEFT JOIN gps.ee_acceleration_limited a' , ...
                    ' ON (s.device_info_serial = a.device_info_serial AND s.date_time = a.date_time)' , ...
                ' LEFT JOIN gps.ee_tracker_limited d ', ...
                    ' ON s.device_info_serial = d.device_info_serial ', ...
                ' LEFT JOIN gps.get_uvagps_track_speed (',num2str(device_numbers(i)), ',' ,start, ',' ,stop, ' ) t ',...
                    ' ON s.device_info_serial = t.device_info_serial and s.date_time = t.date_time ',...
                ' WHERE ( s.n_obs_day >=  ', num2str(minobs),')', ...
                ' ORDER BY s.date_time, a.index' );
        
            % For some odd reason the query only works if converted to char
            sql = char(sql);

            % Execute query and fetch results
            curs = exec(conn, sql);
            disp(curs.Message);
            curs = fetch(curs, batchsize);
            disp(curs.Message);
            
            k = 1;
            while ~strcmp(curs.Data, 'No Data')
                % Storpath needs to exist to store the files
                if exist([storpath,'/', num2str(device_numbers(i))], 'dir') ~= 7
                    mkdir([storpath,'/', num2str(device_numbers(i))]);
                end
                
                data = table2struct(curs.Data, 'ToScalar', true);
                
                % Convert date format for filename:
                % e.g. 2016-08-01 06:00:00 to 20160801
                thisday = dates(j);
                thisday = replace(thisday, '-', ''); % remove dashes
                datechar = char(thisday);
                datechar = datechar(1:8);
                
                % Store file
                filename = sprintf('%s/%dd%sb%d', storpath, device_numbers(i), datechar, k);
                save(filename, 'data');
                
                % Show progress in console
                fprintf('Batch %d of tracker %d on date %s till %s stored \n', k, device_numbers(i), dates(j), dates(j+1));
                
                % Stored, so data can be emptied
                data = [];
                
                % Store fetched device
                devices = [devices device_numbers(i)];
                
                % And fetch a new batch of records
                curs = fetch(curs, batchsize);
                
                k = k + 1;
            end
            
            k = 1;
        end
    end
    
    % All queries done, so we can close the db connection
    close(conn);
    
    % Returns selected devices if successful
    used_devices = unique(devices);
    
    profile viewer;
end

function [datestrings] = generateDateStrings(starttime, endtime, format)
%generateDateStrings Generates date strings for querying data by day
%   Input:
%   1. A starttime in a certain format
%   2. An endtime in a certain format
%   3. The format of starttime, endtime and the query output. If no format
%      is given, the following is used: 'yyyy-mm-dd HH:MM:SS'
%
%   Output:
%   A list of strings iterating throughout the entire datetime range in
%   neat blocks, as follows:
%   2016-08-01 06:00:00
%   2016-08-02 00:00:00
%   2016-08-03 00:00:00
%   2016-08-03 12:00:00

    if nargin < 3
        format = 'yyyy-mm-dd HH:MM:SS';
    end

    startday = datenum(starttime, format);
    startday = floor(startday + 1);
    endday   = datenum(endtime, format);
    endday   = floor(endday);

    datenums = startday:1:endday;
    datestrings = datestr(datenums, format);

    datestrings = string([starttime; datestrings; endtime]);
end

