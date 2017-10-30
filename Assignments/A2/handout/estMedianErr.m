function [F1_noises, median_err] = estMedianErr(pLeft, pRight, nPts,...
                                    xbins, ybins, F0, sigmas, NUM_RESCALE)
    %% Image position noises (Gaussian random noises)

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

    F1_noises = {};
    num_iter = 100;

    for i=1:length(sigmas)
        noise_est = struct();
        noise_est.sigma = sigmas(i);
        for r=1:num_iter
            pLeft_noise = pLeft + [sigmas(i) .* randn(3, nPts)];
            pLeft_noise(3,:) = 1;
            pRight_noise = pRight + sigmas(i) .* randn(3, nPts);
            pRight_noise(3,:) = 1;
            [noise_est.F{r} Sa Sf] = linEstF(pLeft_noise, pRight_noise, NUM_RESCALE);
            noise_est.perpErr(r) = est_perpError(xbins, ybins, ...
                                                F0, noise_est.F{r});
        end
        F1_noises{i} = noise_est;
    end

    median_err = [];
    for i=1:length(sigmas)
       median_err(i) = median(F1_noises{i}.perpErr);
    end


end