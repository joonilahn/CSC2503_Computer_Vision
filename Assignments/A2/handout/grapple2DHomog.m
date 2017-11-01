%% File: grapple2DHomog
%% Uses RANSAC to estimate H matrix from corresponding points.
%% Joonil Ahn

clear
close all
FALSE = 1 == 0;
TRUE = ~FALSE;
global matlabVisRoot

%% We need to ensure the path is set for the iseToolbox.
if length(matlabVisRoot)==0
  dir = pwd;
  cd ../../../Matlab/matlab   %% CHANGE THIS
  startup;
  cd(dir);
end

reconRoot = '.';  %% CHANGE THIS
addpath([reconRoot '/data/wadham']);
addpath([reconRoot '/utils']);


% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

nTrial = 10;  % Number of ransac trials to try.

%% Wadham left image: use  wadham/001.jpg
imPath = 'data/wadham/'; fnameLeft = '001'; 
im = imread([imPath fnameLeft],'jpg');
im = double(im(:,:,2));
imDwn = blurDn(im,1)/2;
imLeft = imDwn;

%% Read correspondence data
load data/wadham/corrPnts5
%% Wadham right image data/wadham002-5.jpg use for corrPnts2-5 respectively
fnameRight = '005';
im = imread([imPath fnameRight],'jpg');
im = double(im(:,:,2));
imDwn = blurDn(im,1)/2;
imRight = imDwn;

clear imPts;
imPts = cat(3, im_pos1', im_pos2');
nPts = size(imPts,2);
if size(imPts,1) == 2
  imPts = cat(1, imPts, ones(1,nPts,2));
end

%% RANSAC for H
seeds = {};
sigma = 2.0; rho = 2;
for kTrial = 1: nTrial
  %% Test out H matrix on a random sample of 4 points
  idTest = randperm(nPts);
  nTest = min(4, nPts);
  idTest = idTest(1:nTest);

  %% Solve for H matrix on the random sample
  [H Sa] = linEstH(imPts(:,idTest,1), imPts(:,idTest,2),1);
  
  %% Compute Euclidean error
  eucErrL = zeros(1,nPts);
  for k = 1:nPts
    left_H = H * imPts(:,k,2);
    left_H = left_H ./ left_H(3);
    eucErrL(k) = norm(left_H - imPts(:,k,1));
  end
  
  %% Detect inliers
  idInlier = abs(eucErrL) < rho*sigma;
  
  %% Count inliers
  nInlier = sum(idInlier);
  if nInlier > 20
    %% Store sets of sampled points with at least 20 inliers
    seed = struct;
    seed.id = idTest;
    seed.idInlier = idInlier;
    seed.nInlier = nInlier;
    seed.H = H;
    
    kSeed = length(seeds)+1;
    seeds{kSeed} = seed;
  end
end 
%% Done RANSAC trials

%% Extract best solution
nInliers = zeros(1, length(seeds));
for ks = 1:length(seeds)
  nInliers(ks) = seeds{ks}.nInlier;
end 
[nM ks] = max(nInliers);
nInliers(ks);

%%  Refine estimate of H using all inliers.
H = seeds{ks}.H;
idInlier = seeds{ks}.idInlier;

idInlierOld = idInlier;
sum(idInlier);
%% Do at most 10 iterations attempting to entrain as many points as possible.
for kIt = 1: 10
  %% Fit F using all current inliers
  [H Sa] = linEstH(imPts(:,idInlier,1), imPts(:,idInlier,2),1);
  
  %% Compute Euclidean error to the transformed points by H
  eucErrL = zeros(1,nPts);
  for k = 1:nPts
    left_H = H * imPts(:,k,2);
    left_H = left_H ./ left_H(3);
    eucErrL(k) = norm(left_H - imPts(:,k,1));
  end
  idInlier = abs(eucErrL) < rho*sigma;
  nInlier = sum(idInlier);
  
  %% If we have the same set of inliers as the previous iteration then stop.
  if all(idInlier == idInlierOld)
    break;
  end
  idInlierOld = idInlier;
end
  

%% Show left image
hFig = figure(1);
clf;
image(imLeft);
colormap(gray(256));
resizeImageFig(hFig, size(imDwn), 1); hold on;
title('Original Left Image');
set(get(hFig,'CurrentAxes'),'Ydir','reverse');

%% Show right image
hFig = figure(2);
clf;
image(imRight);
colormap(gray(256));
resizeImageFig(hFig, size(imDwn), 1); hold on;
title('Original Right Image');
set(get(hFig,'CurrentAxes'),'Ydir','reverse');


%% Warp Images by H
warped_imLeft = homogWarp(imRight, H);
warped_imRight = homogWarp(imLeft, inv(H));

%% Show warped image from left image by H
hFig = figure(3);
clf;
image(warped_imLeft);
colormap(gray(256));
resizeImageFig(hFig, size(imDwn), 1); hold on;
title('Warped Right Image by H');
set(get(hFig,'CurrentAxes'),'Ydir','reverse');

%% Show warped image from right image by H
hFig = figure(4);
clf;
image(warped_imRight);
colormap(gray(256));
resizeImageFig(hFig, size(imDwn), 1); hold on;
title('Warped Left Image by Inverse of H');
set(get(hFig,'CurrentAxes'),'Ydir','reverse');

%% Compute Euclidean distance to warped points in left images.
eucErrL = [];
for k = 1:nPts
  left_H = H * imPts(:,k,2);
  left_H = left_H ./ left_H(3);
  eucErrL(k) = norm(left_H - imPts(:,k,1));
end


%% Plot a histogram of the Euclidean distances
err = eucErrL;
err = min(err, 30);
err = max(err, -30);
[n b] = histo(err, 64);
figure(5); clf;
plot(b,n);
title('Distance from warped image to the original image');
xlabel('Error in pixels');
ylabel('Frequency');

%% Count inliers
idL = abs(eucErrL)< rho*sigma;
% idR = abs(perpErrR) < rho*sigma;
idInlier = idL; %& idR;
sum(idInlier);
sum(idInlier)/nPts;

%% save 'Fcorr' F Sa Sf idInlier nInliers