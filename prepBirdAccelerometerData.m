function [tracks] = prepBirdAccelerometerData(db_user, db_password, device_numbers, starttime, stoptime, minobs, shapefile, storpath)
%prepBirdAccelerometerData Queries accelerometer data
%   Inputs:
%   1. Database username
%   2. Database password
%   3. A vector containing tracker device IDs
%   4. A starttime string in the format 'yyyy-mm-dd hh:mm:ss'
%   5. An endtime string in the format 'yyyy-mm-dd hh:mm:ss'
%   6. A minimum number of observations per day the query will select
%   7. Path to a shapefile delineating the research area bounds

%   Example usage
%   result = prepBirdData(db_user, db_password, [5390, 5419], '2016-08-01 00:00:00', '2016-08-03 00:00:00', 3600, 'data/ResearchAreaNH.shp')
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
    batchsize = 100000; % change also in line below:
    setdbprefs('FetchBatchSize', '100000'); % HAS to be written in quotes for some reason
    setdbprefs('FetchInBatches', 'yes');

    % Turn values into sql-ready strings:
    start = strcat('''', starttime, '''');
    stop = strcat('''', stoptime, '''');
    
    % Determine query boundaries based on the research area shapefile
    shape = shapeinfo(shapefile);
    % Boundingbox contains a 2x2 matrix with minimum and maximum lat and
    % lon values, as follows:
    %         Col 1     Col 2
    % Row 1   min(lon)  min(lat)
    % Row 2   max(lon)  max(lat)
    boundingbox = shape.BoundingBox;
    
    % Setup connection with database
    % Beware: In order to connect, you have to be on the UvA network or
    % using a VPN!
    javaaddpath('drivers/postgresql-42.0.0.jre7.jar'); % required Driver
    conn = database('eecology', db_user, db_password, ...
    'org.postgresql.Driver', 'jdbc:postgresql://db.e-ecology.sara.nl:5432/eecology?ssl=true&sslfactory=org.postgresql.ssl.NonValidatingFactory&');
    setdbprefs('DataReturnFormat', 'table');

    % Query
    % For a description of the table columns, see:
    % https://public.e-ecology.sara.nl/wiki/index.php/DB_Views_2015 and
    % http://nlesc.github.io/eEcology-Classification/#Classification
    tracks = [];
    data   = [];
    
    for i = 1:numel(device_numbers)
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
        disp(curs.Message)

        % Once for a small data sample
        % @TODO: Make code below more elegant and efficient
        curs = fetch(curs, batchsize);
        data = vertcat(data, curs.Data);
        fprintf('First batch of tracker: %d \n', device_numbers(i));
        disp(curs.Message)

        % Now batch the rest
        while ~strcmp(curs.Data, 'No Data')
            curs = fetch(curs, batchsize);
            if ~strcmp(curs.Data, 'No Data')
                data = vertcat(data, curs.Data);
                fprintf('Another batch of: %d records\n', batchsize);
            else
                fprintf('Finished downloading data of tracker: %d \n', device_numbers(i));
            end
        end
        disp(curs.Message);
        
        % Append to tracks and convert data to struct for the classifier
        tracks = vertcat(tracks, data);
        data = table2struct(data, 'ToScalar', true);
        
        % And store it in a .mat file for the classifier
        if exist(storpath, 'dir') ~= 7
            mkdir(storpath);
        end
        
        strbegin = char(data.date_time(1));
        strend   = char(data.date_time(end));
        filename = sprintf('%s/class%d-%still%sminobs%dunannotated.mat', ...
            storpath, device_numbers(i), strbegin(1:10), strend(1:10), minobs);
        save(filename, 'data');
        
        % File saved, now empty data for a next loop iteration
        data = [];
    end
    
    % Close database connection
    close(conn);
    
    profile viewer;
end

