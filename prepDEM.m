function [output_args] = prepDEM(overview, extent, storpath, dimensions, deletefiles)
%prepDEM Downloads and converts DEM data
%   Inputs:
%   1. A shapefile containing an overview of all AHN tiles (ahn_units.shp)
%   2. A shapefile containing the extent of the research area, which will
%      be used to download the corresponding AHN tiles
%   3. A path relative from pwd to store the files. Make sure it exists!
%   4. A 1x2 vector in the form of [width height] where both are in pixels.
%      The raster will be resized to this size.
%   5. Optional: Should unused files be kept? If yes, set to 0. Unused 
%      files will be deleted by default
%
%   Returns 1 if download and conversion is succesful, 0 if not.
%
%   Example use: prepareDEM('data/ahn_units_wgs84.shp','data/ResearchAreaNH2.shp','data/',[2000 2500])
%   
%   Note: Make sure shapefiles are projected in WGS84 and are all stored in
%         the same folder.
%
%   Note: This script may attempt to download non-existing ahn tiles. These
%   will appear as empty .html files in dem-folder. This is because empty
%   tiles (no elevation data available) are not offered. However, these
%   tiles are not removed from the ahn_units shapefile, so the script still
%   attempts to download them.
%   
%   Note: If GDAL commands cannot be found, but you are sure it is
%   installed, paste the value of your PATH variable in the following line
%   of code and unncomment it. Unfortunately MATLAB does not use your
%   system's default PATH.
%   
%   setenv('PATH', ['YOUR_PATH_HERE:', getenv('PATH')])
    
    profile on
    
    if nargin <= 4
        deletefiles = 1; % Default to deleting unused files
    end
    
    % ---------------------------------------------------------------------
    % Clip overview with extent, unless a clipped file already exists in
    % the same location.
    % ---------------------------------------------------------------------
    [pathstr, name, ~] = fileparts(overview);
    clipped = [pathstr,'/',name,'.clipped.shp'];
    
    if exist(clipped,'file') == 0
        % Construct command with following format:
        % ogr2ogr -clipsrc clipfileby.shp output.shp original.shp
        clip = 'ogr2ogr -clipsrc %s %s %s';
        cmd_clip = sprintf(clip, extent, clipped, overview);
        
        [status, cmdout] = system(cmd_clip);
        if status == 1
            error(cmdout)
        end
    end
    
    % ---------------------------------------------------------------------
    % Download AHN files of which units are contained within the clipped
    % shapefile.
    % ---------------------------------------------------------------------
    selection = shaperead(clipped);
    units = string({selection.UNIT});
    nr_units = numel(units);
    
    base_url = 'http://geodata.nationaalgeoregister.nl/ahn2/extract/ahn2_05m_ruw/r%s.tif.zip';
    base_filename = '%sdem/ahn2_05m_ruw_%s.tif.zip';
    dempath = [storpath, 'dem/'];
    
    if exist(dempath, 'dir') ~= 7
        mkdir(dempath); 
    end
    
    parfor i = 1:nr_units
        url = sprintf(base_url,units(i));
        filename = sprintf(base_filename, storpath, units(i));
        
        % First check if a file has already been created before
        file = sprintf('%sr%s.wgs84.tif',dempath,units(i));
        if exist(file,'file') == 2
            fprintf('File %s already exists. Continuing...\n',file);
            continue;
        end
        
        % Download file
        try
            zipped = websave(filename,url);
        catch E
            fprintf('Error requesting: %s\n',url)
            continue;
        end
        
        % Unzipping
        unzipped = unzip(zipped, dempath);
        if deletefiles ~= 0
            delete(zipped);
        end
        
        % Reproject and resize
        % Construct command with following format:
        % gdalwarp pathtosource.tif pathtoreprojected.tif -ts width height ...
        base_reproject = 'gdalwarp "%s" "%sr%s.wgs84.tif" -ts %d %d -r bilinear -t_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"';
        cmd = sprintf(base_reproject, char(unzipped), dempath, units(i), dimensions(1), dimensions(2));
        
        [status,cmdout] = system(cmd);
        if status == 1
            disp(cmdout)
        end
        
        if deletefiles ~= 0
            delete(char(unzipped))
        end
        
        progress = 'Downloaded and reprojected tile %d out of %d. \n';
        fprintf(progress, i, nr_units);
        
    end
    
    profile viewer
    
    output_args = 1;
end