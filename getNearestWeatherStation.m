function [wmoid, distance] = getNearestWeatherStation(stns, lat, lon)
%getNearestWeatherStation Returns the WMO ID of the weather station closest
%   to the point defined by a lat and lon coordinate
%   Input:
%   1. Table with stations information containing a lat, lon and wmo value
%   2. Latitude of point
%   3. Longitude of point
%
%   Output:
%   1. A station ID (wmoid)
%   2. The distance in km from station to the location defined by lat and
%      lon

    n = size(stns, 1);

    distances = zeros(n, 1);

    for i = 1:n
        distances(i) = distWB([stns(i,:).lat stns(i,:).lon], [lat lon]);
    end

    [~, idx] = min(distances);

    wmoid = stns(idx,:).wmo;
    distance = distances(idx);

end

function [km] = distWB(Loc1, Loc2)
    % Loc1 = [Lat1 Lon1]; GPS location A
    % Loc2 = [Lat2 Lon2]; GPS location B
    % example of use to calculate the distance betwee two points in main: 
    % Latitude and Longitude are always presented in degrees.
    % Distance is calculated in [km]
    % Loc1 = [Lat1 Lon1]; GPS location A
    % Loc2 = [Lat2 Lon2]; GPS location B
    % Distance=distWB([52.3 4.5],[52.25 4.58]); or 
    % Distance(i)=distWB([Lat(i) Lon(i)],[Lat(i-1) Lon(i-1)]);
    % W.Bouten, 20160124

    l1=deg2rad(Loc1);
    l2=deg2rad(Loc2);
    lat(1:2)=[l1(1) l2(1)];
    lon(1:2)=[l1(2) l2(2)];

    %% Begin calculation
    R = 6371;                                    % Earth's radius in km
    delta_lat = lat(2) - lat(1);                 % difference in latitude
    delta_lon = lon(2) - lon(1);                 % difference in longitude
    a = sin(delta_lat/2)^2 + cos(lat(1))*cos(lat(2))*sin(delta_lon/2)^2;
    b = 2 * atan2(sqrt(a), sqrt(1-a));
    km = R * b ;                                 % distance in km
end
