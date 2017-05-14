function [deminfo] = loadDEMinfo(storpath, fileformat, nanthreshold)
%loadDEMinfo Returns a table with DEM bounding boxes and summary statistics
%   Inputs:
%   1. A path to the folder containing all the DEM files
%   2. A format to select the right files (e.g. .wgs84.tif or .tif)
%   3. (optional) A threshold value for the raster cells that need to be 
%      converted to NaN. All the values below this value will be set to NaN
%      in order to calculate the accurate summary statistics.
%
%   Output: A table with the following columns:
%   1. The filename
%   2. The min and max X coordinates (longitude)
%   3. The min and max Y coordinates (latitude)
%   4. The data mean, min, max and standard deviation

    files = dir(fullfile(storpath, ['*',fileformat]));

    num_files = size({files.name},2);

    % Pre-alloc
    XMin = zeros(num_files,1);
    XMax = XMin;
    YMin = XMin;
    YMax = XMin;
    geodataMax = XMin;
    geodataMin = XMin;
    geodataMean = XMin;
    geodataStd = XMin;
    
    if nargin > 2
        changeNaN = 1;
    end

    parfor i = 1:num_files

        filename = [storpath,files(i).name];

        info = geotiffinfo(filename);

        % Store tile boundaries
        XMin(i) = info.BoundingBox(1,1);
        XMax(i) = info.BoundingBox(2,1);
        YMin(i) = info.BoundingBox(1,2);
        YMax(i) = info.BoundingBox(2,2);

        % Read file to store summary statistics
        geodata = geotiffread(filename);

        % Set values below threshold to NaN
        if changeNaN == 1
            geodata(geodata <= nanthreshold) = NaN;
        end

        % Store summary statistics
        geodataMax(i) = max(geodata(:));
        geodataMin(i) = min(geodata(:));
        geodataMean(i) = nanmean(geodata(:));
        geodataStd(i) = nanstd(geodata(:));

    end

    filenames = {files.name}';
    deminfo = table(filenames, XMin, XMax, YMin, YMax, geodataMean, geodataMin, ...
                    geodataMax, geodataStd, 'RowNames', filenames);

end

