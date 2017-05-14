clear;
close all;
clc;

profile on

storpath = 'data/dem/';

files = dir(fullfile(storpath, '*.wgs84.tif'));

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
    geodata(geodata < -1e10) = NaN;
    
    % Store summary statistics
    geodataMax(i) = max(geodata(:));
    geodataMin(i) = min(geodata(:));
    geodataMean(i) = nanmean(geodata(:));
    geodataStd(i) = nanstd(geodata(:));
    
end

filenames = {files.name}';
deminfo = table(filenames, XMin, XMax, YMin, YMax, geodataMean, geodataMin, ...
                geodataMax, geodataStd, 'RowNames', filenames);

profile viewer