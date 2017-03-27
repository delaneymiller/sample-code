%% Quantifying Particle Adherence to Acidic PEG Fibers (Multiple Images)
% 29 July 2016
% Author: Delaney Miller
% Contributors: Hyoungsoo Kim, Janine Nunes
% Complex Fluids Group, Princeton University

%% Image Processing Procedure
% 1) Initialize "for" loop to read through all the .tif files in a folder
% 2) Create binary mask
% 3) Find area of binary mask (in pixels)
% 4) Apply binary mask to fluorescent image
% 5) Binarize the fluorescent image and save to folder
% 6) Measure and count the objects

%% Preliminary Steps

% Initialize the memory
close all;
clear all; clc;
close all;
clear all;

% Set high precision
format LONGENG
               
%% 1) Initialize "for" loop to read through all the .tif files in a folder

pathname = '/Users/DMiller/Documents/MATLAB/Fiber Data/25% acid';
cd(pathname);
grayFiles = dir('*ch01.tif');
fluoFiles = dir('*ch00.tif');
numImages = length(grayFiles(not([grayFiles.isdir])));

% Make a new folder for the post-processed images
f2 = fullfile(pathname,'25% processed images');
if (exist(f2) == 0)
    mkdir(f2);  % create folder
end

% Create an array to store the area of the fiber in each image
fiberAreas = zeros(numImages,1);

% Create an array to store the number of objects in each image
objects = zeros(numImages,1);
objectAreas = [];

% Initialize "for" loop
for k = 1 : length(grayFiles)
    im1 = imread(grayFiles(k).name);
    im1 = im1(:,:,1);
    imf1 = imread(fluoFiles(k).name);
    imf2 = imf1(:,:,1);

    % Enhance contrast
    im2 = imadjust(im1); 
    
    % Find edges of objects in the image (binarize)
    im3 = edge(im2,'sobel','nothinning');
    im4a = imdilate(im3,strel('disk',2));
    im4b = imerode(im4a,strel('disk',2));
    im4 = imfill(im4b,'holes');

    % Add borders to top and bottom
    im4(1,:) = 1;
    im4(end,:) = 1;

    % Close all object holes
    im4 = imfill(im4,'holes');

    % Find area of largest object
    stats1 = regionprops(im4);
    areas1 = [stats1.Area];
    max1 = max(areas1);

    % Remove all objects from image except the largest
    im5 = bwareafilt(im4,[max1 - 1, max1 + 1]);

    % Get rid of borders
    im5(1,:) = 0;
    im5(end,:) = 0;

	% Erase additional objects
    stats2 = regionprops(im5);
    areas2 = [stats2.Area];
    max2 = max(areas2);
    im6 = bwareafilt(im5,[max2 - 1, max2 + 1]);

    % Dilate and erode fiber to smooth out
    im7 = imdilate(im6,strel('disk',2));
    im8 = imerode(im7,strel('disk',2));

    %% 3) Find area of binary mask
    fiberAreas(k,:) = bwarea(im8);

    %% 4) Apply binary mask to fluorescent image
    for i = 1 : size(im8,2)
        for j = 1 : size(im8,1)
            if im8(i,j) == 0
                imf2(i,j) =0;
            end
        end
    end

    %% 5) Binarize the fluorescent image and save to folder

    % Apply Gaussian smoothing filter to detect particle signal   
    imf3 = imgaussfilt(imf2,0.3);
    
    % Binarize the image
    imf4 = im2bw(imf3,0.4);
 
    % Save the mask and image to the processed images folder
    maskName = sprintf('mask%d',k);
    resultName = sprintf('result%d', k);
    imwrite(im8,[f2,'/',maskName,'.tif'],'Compression','none');
    imwrite(imf4,[f2,'/',resultName,'.tif'],'Compression','none');

    %% 6) Measure and count the objects
    % Measure the areas of the objects in the image
    stats  = regionprops(imf4, 'Area');
    objectAreas = cat(1, objectAreas, stats.Area);
    
    % Count how many objects are in the image and add to array
    objects(k,:) = length(stats);
end

% total objects
totObjects = sum(objects)

% Find total fiber area (sum from all images)
totFiberArea = sum(fiberAreas)

% Find total object area (sum from all images)
totObjectArea = sum(objectAreas)

% Find standard deviation of objects
std = std(objects)

%% Make a histogram to determine actual number of particles
histo = histogram(objectAreas,'BinMethod','Integers');
title('25% Acid Fibers: Histogram of Object Sizes')
xlabel('Pixels')
ylabel('Number of objects')
[histo1, edges] = histcounts(objectAreas,'BinMethod','integers');

% Find size and index of maximum bin
[maxV, maxI] = max(histo1);
bin = [edges(maxI), edges(maxI + 1)];

% Find average particle size
avgParticle = mean(bin);

% Find total number of particles based on average particle size
totParticle = totObjectArea/avgParticle;

% Normalize this w/r/t the total fiber area
totObjectDensity = totObjectArea / totFiberArea