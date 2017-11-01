function median_err = estMedianErr(pLeft, pRight, nPts,...
                                    xbins, ybins, F0, sigmas, NUM_RESCALE)
    %% Image position noises (Gaussian random noises)

    if nargin < 8
        NUM_RESCALE = 1;
    end

    F1_noises = {};
    num_iter = 100;

    for i=1:length(sigmas)
        noise_est = struct();
        noise_est.sigma = sigmas(i);
        for r=1:num_iter
            pLeft_noise = pLeft + sigmas(i) .* randn(3, nPts);
            pLeft_noise(3,:) = 1;
            pRight_noise = pRight + sigmas(i) .* randn(3, nPts);
            pRight_noise(3,:) = 1;         
            F1 = linEstF(pLeft_noise, pRight_noise, NUM_RESCALE);
            noise_est.perpErr(r) = est_perpError(xbins, ybins, F0, F1);
        end
        F1_noises{i} = noise_est;
    end

    median_err = [];
    for i=1:length(sigmas)
       median_err(i) = median(F1_noises{i}.perpErr);
    end


end