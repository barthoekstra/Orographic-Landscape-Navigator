function [tracks] = prepBirdData(db_user, db_password, device_numbers, starttime, stoptime, shapefile)
%prepBirdData Queries the database and returns corresponding data
%   Inputs:
%   1. Database username
%   2. Database password
%   3. A vector containing tracker device IDs
%   4. A starttime string in the format 'yyyy-mm-dd hh:mm:ss'
%   5. An endtime string in the format 'yyyy-mm-dd hh:mm:ss'
%   6. Path to a shapefile delineating the research area bounds

%   Example usage
%   result = prepBirdData(db_user, db_password, [5390, 5419], '2016-08-01 00:00:00', '2016-08-03 00:00:00', 'data/ResearchAreaNH.shp')
%
%   All tracking data in a table, with device serial, date and time,
%   longitude, latitude, altitude, temperature, gps time to fix, 2d speed,
%   altitude above ground level, positional dilution of precision,
%   satellites used, nr of acceleration measurements, z-axis acceleration
%   variance and the project the tracker belongs to.
%   

    % Batching
    % Although not necessary, the code is written to allow for batching the
    % download. I have not checked to what extent I can push the limits of
    % this code until I encounter memory issues, but the following
    % variables can be changed to do so. More info on the url below:
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
    % https://public.e-ecology.sara.nl/wiki/index.php/DB_Views_2015
    %
    % Query adapted from LBBG_Texel2 code by Willem Bouten and adjusted to
    % only select records within the study area and to include details on
    % the project.
    tracks = [];
    
    for i = 1:numel(device_numbers)
        sql = strcat('SELECT t.device_info_serial, t.date_time, ', ...
                ' date_part(''year''::text, t.date_time) AS year, ', ...
                ' date_part(''month''::text, t.date_time) AS month, ', ...
                ' date_part(''day''::text, t.date_time) AS day, ', ...
                ' date_part(''hour''::text, t.date_time) AS hour, ', ...
                ' date_part(''minute''::text, t.date_time) AS minute, ', ...
                ' date_part(''second''::text, t.date_time) AS second, ', ...
                ' COUNT(*) OVER (PARTITION BY DATE(t.date_time)) AS n_obs_day, ', ...
                ' t.longitude, t.latitude, t.altitude, t.temperature, t.gps_fixtime, ', ...
                ' t.speed_2d, t.altitude_agl, ', ...
                ' t.positiondop, t.satellites_used, ', ...
                ' p.key_name AS project', ...
            ' FROM gps.ee_tracking_speed_limited t' , ...
            ' LEFT JOIN gps.ee_tracker_ownership_limited p', ...
            ' ON (t.device_info_serial = p.device_info_serial)', ...
            ' WHERE ( t.date_time >= ', start,  ')', ...
                ' AND ( t.date_time <= ', stop,    ')' , ...
                ' AND ( t.device_info_serial = ', num2str(device_numbers(i)), ')', ...
                ' AND ( t.longitude >= ', num2str(boundingbox(1,1)), ')', ...
                ' AND ( t.longitude <= ', num2str(boundingbox(2,1)), ')', ...
                ' AND ( t.latitude >= ', num2str(boundingbox(1,2)), ')', ...
                ' AND ( t.latitude <= ', num2str(boundingbox(2,2)), ')', ...
                ' AND t.userflag = 0', ...
            ' GROUP BY p.key_name,  t.device_info_serial, t.date_time, t.longitude, ' , ...
                ' t.speed_2d, t.altitude_agl,' , ...
                ' t.latitude, t.altitude, t.temperature, t.gps_fixtime, t.positiondop, t.satellites_used, t.location ' , ...
            ' ORDER BY t.date_time' );

        % For some odd reason the query only works if converted to char
        sql = char(sql);

        % Execute query and fetch results
        curs = exec(conn, sql);

        % Once for a small data sample
        % @TODO: Make code below more elegant and efficient
        curs = fetch(curs, batchsize);
        tracks = vertcat(tracks, curs.Data);
        fprintf('First batch of: %d \n', device_numbers(i));
        disp(curs.Message)

        % Now batch the rest
        k = 1;
        while ~strcmp(curs.Data, 'No Data')
            curs = fetch(curs, batchsize);
            if ~strcmp(curs.Data, 'No Data')
                tracks = vertcat(tracks, curs.Data);
                fprintf('Batch %d of tracker %d \n', k, device_numbers(i));
                k = k + 1;
            else
                fprintf('Finished downloading data of: %d \n', device_numbers(i));
                k = 1;
            end
        end
        disp(curs.Message);
    end
    
    % Close database connection
    close(conn);
    
end

