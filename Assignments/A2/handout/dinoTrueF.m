% File: dinoTrueF.m
% Script for estimation of True F (F0) and compare with F1

clear
close all
FALSE = 1 == 0;
TRUE = ~FALSE;

global matlabVisRoot

%% We need to ensure the path is set for the iseToolbox.
if length(matlabVisRoot)==0
  dir = pwd;
  % cd ~jepson/pub/matlab   %% CHANGE THIS to your startup directory
  cd /Users/Joonil/Documents/2017 Fall/CSC2503_Foundations of Computer Vision/Matlab/matlab
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
nInliers(ks)

%%  Refine estimate of F using all inliers.
F1 = seeds{ks}.F;
idInlier = seeds{ks}.idInlier;

idInlierOld = idInlier;
sum(idInlier)
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
Tvec = dLeft - dRight;
T = zeros(3,3);
T(1) = 0;
T(2) = Tvec(3);
T(3) = Tvec(2);
T(4) = Tvec(3);
T(5) = 0;
T(6) = -Tvec(1);
T(7) = -Tvec(2);
T(8) = Tvec(1);
T(9) = 0;
E = Rleft* T * Rright';                         % Essential matrix
F0 = inv(MintLeft)' * E * inv(MintRight);     % Ground-truth F matrix


%% A.2. Compare F0 and F1
nPoints = 10;  %% Number of points per axis
nCol = 8;   %% Number of different colours to use.
col = hsv(nCol);  %% Colour map.

figure(1);
xbins = linspace(-150, 150, nPoints);
ybins = linspace(-100, 100, nPoints);
% xbins = linspace(0, 100, nPoints);
% ybins = linspace(0, 100, nPoints);
[X, Y] = meshgrid(xbins, ybins);
plot(X, Y, '.b');
axis([-150 150 -100 100]); axis xy;
% axis([0 100 0 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Regularly spaced grid of points');

hfig = figure(2);
clf; hold on;
title('Epipolar Lines for the grid of points');
axis([-150 150 -100 100]);
% axis([0 100 0 100]); 
% axis xy; axis equal;
perpErr = [];
for i = 1:nPoints
    for j = 1:nPoints
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

sigmas = [0.001 0.01 0.1 1.0];
F1s_noises = {};

for i=1:length(sigmas)
    noise_est = struct();
    noise_est.sigma = sigmas(i);
    for r=1:100
        pLeft_noise = pLeft + randn(3, nPts) * sigmas(i);
        pLeft_noise(3,:) = 1;
        pRight_noise = pRight + randn(3, nPts) * sigmas(i);
        pRight_noise(3,:) = 1;
        [noise_est.F{r} Sa Sf] = linEstF(pLeft_noise, pRight_noise, 1);
        noise_est.perpErr(r) = est_perpError(xbins, ybins, ...
                                             F0, noise_est.F{r}, nPoints);
    end
    F1_noises{i} = noise_est;
end

median_err = [];
for i=1:length(sigmas)
   median_err(i) = median(F1_noises{i}.perpErr);
end
figure(3);
axis([0.0 1.0 0 15]); axis xy;
semilogx(sigmas, median_err);
title('Median error as a function of sigma');
xlabel('log sigma');
ylabel('Median of maximum perpendicular error');

%% A.4. No hartley's normalization
for i=1:length(sigmas)
    noise_est = struct();
    noise_est.sigma = sigmas(i);
    for r=1:100
        pLeft_noise = pLeft + randn(3, nPts) * sigmas(i);
        pLeft_noise(3,:) = 1;
        pRight_noise = pRight + randn(3, nPts) * sigmas(i);
        pRight_noise(3,:) = 1;
        [noise_est.F{r} Sa Sf] = linEstF(pLeft_noise, pRight_noise);
        noise_est.perpErr(r) = est_perpError(xbins, ybins, ...
                                             F0, noise_est.F{r}, nPoints);
    end
    F1_noises2{i} = noise_est;
end

median_err2 = [];
for i=1:length(sigmas)
   median_err2(i) = median(F1_noises2{i}.perpErr);
end
figure(4);
semilogx(sigmas, median_err2);
title('Median error as a function of sigma with no Hartley''s Normalization');
xlabel('log sigma');
ylabel('Median of maximum perpendicular error');

%% A.5. Change the scale of Z component of Dino
sclZ = 0.1;
% Dino left image data
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, Rleft, sclZ);

% Dino right image data
[pRight polys MintRight MextRight] = projectDino(f, dRight, Rright, sclZ);

for i=1:length(sigmas)
    noise_est = struct();
    noise_est.sigma = sigmas(i);
    for r=1:100
        pLeft_noise = pLeft + randn(3, nPts) * sigmas(i);
        pLeft_noise(3,:) = 1;
        pRight_noise = pRight + randn(3, nPts) * sigmas(i);
        pRight_noise(3,:) = 1;
        [noise_est.F{r} Sa Sf] = linEstF(pLeft_noise, pRight_noise, 1);
        noise_est.perpErr(r) = est_perpError(xbins, ybins, ...
                                             F0, noise_est.F{r}, nPoints);
    end
    F1_noises3{i} = noise_est;
end

median_err3 = [];
for i=1:length(sigmas)
   median_err3(i) = median(F1_noises3{i}.perpErr);
end
figure(5);
semilogx(sigmas, median_err3);
title('Median error as a function of sigma with sclZ=0.1');
xlabel('log sigma');
ylabel('Median of maximum perpendicular error');
