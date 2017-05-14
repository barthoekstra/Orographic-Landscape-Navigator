%% Creating a Half-Resolution Georeferenced Image
%
% This example shows how to create a half-resolution version of a
% georeferenced TIFF image, using referencing objects and Image Processing
% Toolbox(TM) functions |ind2gray| and |imresize|.

% Copyright 1996-2014 The MathWorks, Inc.

%% Step 1: Import a Georeferenced TIFF Image
% Read an indexed-color TIFF image and convert it to grayscale.
% The size of the image is 2000-by-2000.
[X, cmap, R] = geotiffread('data/ahn/r14bn2.tif');
% I_orig = ind2gray(X, cmap);

%%
% Read the corresponding worldfile.  Each image pixel covers a
% one-meter square on the map.  
R_orig = worldfileread('data/ahn/r14bn2.tif','planar',size(X));

%%
% Choose a convenient format for displaying the result.
currentFormat = get(0,'format');
format short g
R_orig

%% Step 2: Resize the Image to Half Its Original Size
% Halve the resolution, creating a smaller (1000-by-1000) image.
I_half = imresize(I_orig, size(I_orig)/2, 'bicubic');

%% Step 3: Construct a Referencing Object for the Resized Image
% The resized image has the same limits as the original, just a smaller
% size, so copy the referencing object and reset its RasterSize property.

R_half = R_orig;
R_half.RasterSize = size(I_half)

%% Step 4: Visualize the Results
% Display each image in map coordinates, and mark a reference point with a
% red + in both figures.

xlimits = [208000 208250];
ylimits = [911800 911950];

x = 208202.21;
y = 911862.70;

figure
mapshow(I_orig,R_orig)
hold on
plot(x,y,'r+')
xlim(xlimits)
ylim(ylimits)
ax = gca;
ax.TickDir = 'out';

%%
figure
mapshow(I_half,R_half)
hold on
plot(x,y,'r+')
xlim(xlimits)
ylim(ylimits)
ax = gca;
ax.TickDir = 'out';

%%
% Graphically, they coincide, even though the same map location
% corresponds to two different locations in intrinsic coordinates.
[xIntrinsic1, yIntrinsic1] = worldToIntrinsic(R_orig, x, y)

%%
[xIntrinsic2, yIntrinsic2] = worldToIntrinsic(R_half, x, y)

format(currentFormat);

%% Credits
% concord_ortho_w.tif, concord_ortho_w.tfw - derived from orthophoto
% tiles from:
%
%    Office of Geographic and Environmental Information (MassGIS),
%    Commonwealth of Massachusetts  Executive Office of Environmental Affairs
%    http://www.state.ma.us/mgis
%    
%    For more information, run: 
%    
%    >> type concord_ortho.txt
%
