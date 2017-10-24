function filt = catmullRom(subSample)
% Return the Catmull-Rom filter kernel, subsampled by subSample.
% filt will be 2*(2*subSample-1) + 1 taps long.

  x = (1:subSample)/subSample;
 
  filt = 1 + 0.5 * (-5 * x .*x + 3 * x .* x .* x);
  filt(end) = 0.0;
  filt = [filt 0.5*(-x + 2* x .* x - x.*x.*x)];
  filt = filt(1:(end-1));
  filt = [filt(end:-1:1) 1 filt];
