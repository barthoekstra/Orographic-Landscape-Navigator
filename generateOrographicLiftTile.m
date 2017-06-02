function [ OrographicLift ] = generateOrographicLiftTile(demfile, nanthreshold, cellsize, wspeed, wdir, storname)
%generateOrographicLiftTile Generates orographic lift values for given tile
%   This script quickly generates orographic lift values for a given DEM
%   tile and optionally stores this as a resulting GeoTIFF file
%   
%   Input:
%   1. GeoTIFF DEM file
%   2. NaN threshold (values below will be set to 0)
%   3. Cellsize in relevant units
%   4. Wind speed
%   5. Wind direction
%   6. Output filename

    if nargin > 5
        storefile = 1;
    else
        storefile = 0;
    end

    dem = geotiffread(demfile);
    deminfo = geotiffinfo(demfile);
    dem(dem <= nanthreshold) = 0;

    [r, c] = size(dem);

    % dzdy
    dzdy = zeros(r, c);
    dzdy(1,:) = (dem(1,:) - dem(2,:)) ./ (cellsize);
    dzdy(r,:) = (dem(r-1,:) - dem(r,:)) ./ (cellsize);
    dzdy(2:r-1,:) = (dem(1:r-2,:) - dem(3:r,:)) ./ (cellsize * 2);

    % dzdx
    dzdx = zeros(r, c);
    dzdx(:,1) = (dem(:,1) - dem(:,2)) ./ (cellsize);
    dzdx(:,c) = (dem(:,c-1) - dem(:,c)) ./ (cellsize);
    dzdx(:,2:c-1) = (dem(:,1:c-2) - dem(:,3:c)) ./ (cellsize * 2);

    % Slope
    slope = atand(sqrt(dzdx.^2 + dzdy.^2));

    % Aspect
    aspect = 90 + atan2d(dzdy, dzdx);
    aspect = aspect + (aspect < 0) * 360; % transform to [0 360]

    % Cleanup
    slope(isnan(slope)) = 0;

    Ca = sind(slope) .* cosd(wdir - aspect);

    OrographicLift = wspeed .* Ca;

    if storefile == 1
        geotiffwrite(storname, OrographicLift, deminfo.SpatialRef)
    end

end
