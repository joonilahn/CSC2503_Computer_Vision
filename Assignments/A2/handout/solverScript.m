% File: dinoTrueF.m
% Script for estimation of True F (F0) and compare with F1
% Author: Joonil Ahn

clear
close all
FALSE = 1 == 0;
TRUE = ~FALSE;

global matlabVisRoot

%% We need to ensure the path is set for the iseToolbox.
if length(matlabVisRoot)==0
  dir = pwd;
  % cd ~jepson/pub/matlab   %% CHANGE THIS to your startup directory
  cd ../../../Matlab/matlab
  startup;
  cd(dir);
end

reconRoot = '.';  %% CHANGE THIS to the directory you installed A4 handout/
addpath([reconRoot '/utils']);


% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

nTrial = 10;  %% Number of ransac trials to use


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up cameras
%% The cameras are automatically rotated by projectDino to fixate
%% on the mean of the 3D points.  We do not always want to allow
%% the cameras to move in this fashion (eg. when we change sclZ).
%% So we will compute the rotations for the left and right cameras
%% once and for all, and then use these.
f = 100; % focal length
dLeft = [-50, 0, -150]';  % Center of projection for left camera
dRight = [50, 0, -150]';  % Center of projection for right camera
%% Compute camera rotations to fixate on Dino's center.
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, [], 1.0);
Rleft = MextLeft(:, 1:3);
[pRight polys MintRight MextRight] = projectDino(f, dRight, [], 1.0);
Rright = MextRight(:, 1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate data...
sclZ = 1.0;
%% Dino left image data
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, Rleft, sclZ);

%% Dino right image data
[pRight polys MintRight MextRight] = projectDino(f, dRight, Rright, sclZ);

%% Build correspondence data
clear imPts;
imPts = cat(3, pLeft, pRight);
nPts = size(imPts,2);
if size(imPts,1) == 2
  imPts = cat(1, imPts, ones(1,nPts,2));
end

%% RANSAC for F
seeds = {};
sigma = 2.0; rho = 2;
for kTrial = 1: nTrial
  %% Test out F matrix on a random sample of 8 points
  idTest = randperm(nPts);
  nTest = min(8, nPts);
  idTest = idTest(1:nTest);

  %% Solve for F matrix on the random sample
  [F1 Sa Sf] = linEstF(imPts(:,idTest,1), imPts(:,idTest,2),1);
  
  %% Compute perpendicular error of all points to epipolar lines
  perpErrL = zeros(1,nPts);
  for k = 1:nPts
    lk = imPts(:,k,2)' * F1';
    perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
  end
  
  %% Detect inliers
  idInlier = abs(perpErrL) < rho*sigma;
  
  %% Count inliers
  nInlier = sum(idInlier);
  if nInlier > 20
    %% Store sets of sampled points with at least 20 inliers
    seed = struct;
    seed.id = idTest;
    seed.idInlier = idInlier;
    seed.nInlier = nInlier;
    seed.F = F1;
    
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

%%  Refine estimate of F using all inliers.
F1 = seeds{ks}.F;
idInlier = seeds{ks}.idInlier;

idInlierOld = idInlier;

%% Do at most 10 iterations attempting to entrain as many points as possible.
for kIt = 1: 10
  %% Fit F using all current inliers
  [F1 Sa Sf] = linEstF(imPts(:,idInlier,1), imPts(:,idInlier,2),1);
  
  %% Compute perpendicular error to epipolar lines
  perpErrL = zeros(1,nPts);
  for k = 1:nPts
    lk = imPts(:,k,2)' * F1';
    perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
  end
  idInlier = abs(perpErrL) < rho*sigma;
  nInlier = sum(idInlier);
  
  %% If we have the same set of inliers as the previous iteration then stop.
  if all(idInlier == idInlierOld)
    break;
  end
  idInlierOld = idInlier;
end
  
%% A.1. Compute true fundamental matrix F0
F0 = estTrueF(dLeft, dRight, Rleft, Rright, MintLeft, MintRight);

%% A.2. Compare F0 and F1
xinterval = 31;  %% Number of points per x-axis
yinterval = 21;  %% Number of points per y-axis
nCol = 8;   %% Number of different colours to use.
col = hsv(nCol);  %% Colour map.

figure(1);
xbins = linspace(-150, 150, xinterval);
ybins = linspace(-100, 100, yinterval);
[X, Y] = meshgrid(xbins, ybins);
plot(X, Y, '.b');
axis([-150 150 -100 100]); axis xy;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Regularly spaced grid of points');

figure(2);
title('Epipolar Lines for the grid of points');
% axis([-150 150 -100 100]);
perpErr = [];
for i = 1:xinterval
    for j = 1:yinterval
      % Plot interest point location corresponding to epipolar line as a "o"
      % in the same colour as the epipolar line.
%       plot(xbins(i), ybins(j), 'o', 'Color', 'b');
      % Plot epipolar line.
      lF0 = [xbins(i); ybins(j); 1]' * F0';
      lF1 = [xbins(i); ybins(j); 1]' * F1';
      ep0 = cropLineInBox(lF0(1:2), lF0(3), cropBox);
      if sum(sum(isnan(ep0))) > 0
          break;
      else
        set(line(ep0(:,1), ep0(:,2)), 'Color',  col(mod(i,nCol)+1,:));
        maxperpErr_ij = max(abs([ep0 ones(2,1)]*lF1') / norm(lF1(1:2)));
        perpErr = [perpErr maxperpErr_ij];
      end
    end
end

max_perpErr = max(perpErr);     % compute max of max perpendicular error
fprintf('The maximum of maximum perpendicular errors is %e\n', max_perpErr);

%% A.3. Image position errors (Gaussian random noises)

% % Random number generator seed:
% seed = round(sum(1000*clock));
% rand('state', seed);
% seed0 = seed;
% % We also need to start randn. Use a seedn generated from seed:
% seedn = round(rand(1,1) * 1.0e+6);
% randn('state', seedn);

sigmas = [0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1.0 2.0];
[F1_noises3, median_err3] = estMedianErr(pLeft, pRight, nPts, xbins, ybins,...
                                                        F0, sigmas);

% %% A.4. No hartley's normalization
% % Estimate median errors
% [F1_noises4, median_err4] = estMedianErr(pLeft, pRight, nPts, xbins, ybins,...
%                                                         F0, sigmas, 0);
% 
% %% A.5. Change the scale of Z component of Dino
% f = 100;
% sclZ = 0.1;
% 
% % Dino left image data with sclZ=0.1
% [pLeft_sclZ polys MintLeft MextLeft] = projectDino(f, dLeft, Rleft, sclZ);
% % Dino right image data with sclZ=0.1
% [pRight_sclZ polys MintRight MextRight] = projectDino(f, dRight, Rright, sclZ);
% 
% % Reestimate the F0
% F0_sclZ = estTrueF(dLeft, dRight, Rleft, Rright, MintLeft, MintRight);                                  
% 
% % Estimate median errors
% [F1_noises5, median_err5] = estMedianErr(pLeft_sclZ, pRight_sclZ, nPts,...
%                                       xbins, ybins, F0_sclZ, sigmas);

%% A.6. Smaller separation between two cameras
% Different locations of centre of projection for the two cameras
sclZ = 1;
dLeft_sm = [0, 0, -150]';
dRight_sm = [0, 0, -150]';

% Compute camera rotations for the different separation
[pLeft_sm polys MintLeft_sm MextLeft_sm] = projectDino(f, dLeft_sm, [], 1);
Rleft_sm = MextLeft_sm(:, 1:3);
[pRight_sm polys MintRight_sm MextRight_sm] = projectDino(f, dRight_sm, [], 1);
Rright_sm = MextRight_sm(:, 1:3);

% Dino left image data
[pLeft_sm polys MintLeft_sm MextLeft_sm] = projectDino(f, dLeft_sm, Rleft_sm, 1);
% Dino right image data
[pRight_sm polys MintRight_sm MextRight_sm] = projectDino(f, dRight_sm, Rright_sm, 1);

% Re-estimate the F0
F0_sm = estTrueF(dLeft_sm, dRight_sm, Rleft_sm, Rright_sm, MintLeft_sm, MintRight_sm);

% Estimate median errors
[F1_noises6, median_err6] = estMedianErr(pLeft_sm, pRight_sm, nPts,...
                                      xbins, ybins, F0_sm, sigmas);

%% Plot figures for A.3 A.4 A.5 A.6         
% Plot figures for A.3
figure(3); hold on;
plot(sigmas, median_err3);
set(gca, 'xscale','log')
% axis tight;
title('Median error as a function of sigma');
xlabel('log sigma');
ylabel('Median of maximum perpendicular error');

% % Plot figures for A.4
% figure(4);
% plot(sigmas, median_err4);
% set(gca, 'xscale','log')
% title('Median error as a function of sigma with no Hartley''s Normalization');
% xlabel('log sigma');
% ylabel('Median of maximum perpendicular error');
% 
% % Plot figure for A.5
% figure(5);
% plot(sigmas, median_err5);
% set(gca, 'xscale','log')
% title('Median error as a function of sigma with sclZ=0.1');
% xlabel('log sigma');
% ylabel('Median of maximum perpendicular error');

% Plot figure for A.6
figure(6);
plot(sigmas, median_err6);
set(gca, 'xscale','log')
title('Median error as a function of sigma with smaller separation');
xlabel('log sigma');
ylabel('Median of maximum perpendicular error');

