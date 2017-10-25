function [F1_noises, median_err] = estMedianErr(pLeft, pRight, nPts,...
                                    xbins, ybins, F0, sigmas, NUM_RESCALE)
    %% A.3. Image position errors (Gaussian random noises)

    % % Random number generator seed:
    % seed = round(sum(1000*clock));
    % rand('state', seed);
    % seed0 = seed;
    % % We also need to start randn. Use a seedn generated from seed:
    % seedn = round(rand(1,1) * 1.0e+6);
    % randn('state', seedn);
    if nargin < 8
        NUM_RESCALE = 1;
    end

    F1s_noises = {};
    pts_axis = size(xbins, 1);
    num_iter = 100;

    for i=1:length(sigmas)
        noise_est = struct();
        noise_est.sigma = sigmas(i);
        for r=1:100
            pLeft_noise = pLeft + randn(3, nPts) * sigmas(i);
            pLeft_noise(3,:) = 1;
            pRight_noise = pRight + randn(3, nPts) * sigmas(i);
            pRight_noise(3,:) = 1;
            [noise_est.F{r} Sa Sf] = linEstF(pLeft_noise, pRight_noise, NUM_RESCALE);
            noise_est.perpErr(r) = est_perpError(xbins, ybins, ...
                                                 F0, noise_est.F{r}, pts_axis);
        end
        F1_noises{i} = noise_est;
    end

    median_err = [];
    for i=1:length(sigmas)
       median_err(i) = median(F1_noises{i}.perpErr);
    end


end