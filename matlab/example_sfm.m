% EXAMPLE_SFM Mini tutorial on SfM from two images
%
% Other m-files required: ./@Matcha/*, ./geo/*, ./math/*, ./3d/*, ./graphix/*
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

close all;
clear all;
clc;

%% Import the good stuff

addpath('./geo');
addpath('./3d');
addpath('./math');
addpath('./graphix');

%% Feature extraction image 1

exif1 = imfinfo("./data/doge_palace/DSC_0059.JPG");
rgb1  = imread("./data/doge_palace/DSC_0059.JPG");

im1 = Matcha(rgb1);
im1.camera.K = intrinsicsMatrix(exif1.DigitalCamera.FocalLength / 0.00605, ...
                                exif1.DigitalCamera.FocalLength / 0.00605, ...
                                size(rgb1, 2), size(rgb1, 1), 0);

im1.params.maxOctave = 4;
im1.params.edgeThreshold  = 0.06;
im1.params.edgeThreshold1 = 0.06;
im1.params.peakThreshold  = 255 * 0.04 / 3.0;

im1.computePyramid();
im1.computeDogExtrema();
im1.filterEdges();
im1.interpolateExtrema();
im1.computeOrientations();
im1.computeFeatures();

%% Feature extraction image 2

exif2 = imfinfo("../data/doge_palace/DSC_0059.JPG");
rgb2  = imread("../data/doge_palace/DSC_0061.JPG");
im2   = Matcha(rgb2);

im2.camera.K = intrinsicsMatrix(exif2.DigitalCamera.FocalLength / 0.00605, ...
                                exif2.DigitalCamera.FocalLength / 0.00605, ...
                                size(rgb2, 2), size(rgb2, 1), 0);

im2.params.maxOctave = 4;
im2.params.edgeThreshold  = 0.06;
im2.params.edgeThreshold1 = 0.06;
im2.params.peakThreshold  = 255 * 0.04 / 3.0;

im2.computePyramid();
im2.computeDogExtrema();
im2.filterEdges();
im2.interpolateExtrema();
im2.computeOrientations();
im2.computeFeatures();

%% Draw keypoints

subplot(1,2,1);
im1.showExtrema();

subplot(1,2,2);
im2.showExtrema();

%% Match nearest neighbors (im1 with im2)

[inliers1, inliers2, F12] = im1.match(im2, 0.75, 1, 2000, 0.02, 50);

%% Eliminate large residuals

residuals = diag([inliers1, ones(size(inliers1,1),1)] * F12 * [inliers2'; ones(1,size(inliers2,1))]);

stddev = std(residuals);
idx = find(residuals(abs(residuals) < 2 * stddev,:));

inliers1 = inliers1(idx,:);
inliers2 = inliers2(idx,:);

F12 = epipolarGeometry(inliers1, inliers2);

%% Show matches
            
showMatches(im1.rgb, inliers1, im2.rgb, inliers2, 'side', 1.0);

%% 3D Reconstruction

% Essential matrix
E12 = im1.camera.K' * F12 * im2.camera.K;

% Extract pose from essential
[R1, t1, R2, t2] = computePose(E12);

% Create projection matrices
P1 = projectionMatrix(im1.camera.K, eye(3), zeros(3,1));
P2 = projectionMatrix(im2.camera.K, R1', t2);

pcloud.pts   = triangulate(inliers1, inliers2, P1, P2);
pcloud.color = im1.at(inliers1);
    
plotCloud(pcloud, 15);
