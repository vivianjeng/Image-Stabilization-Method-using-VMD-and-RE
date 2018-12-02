function [g,x] = gaussian_filter( sigma, sample )

% g = gaussian_filter( sigma )
%
% Create 1D a gaussian filter with a suitable number 
% of filter taps for the specified sigma.
%
% Input
% sigma - standard deviation of the gaussian
% sample - multiple of sigma to sample filter to
%
% Output
% g - the gaussian filter
%
% Thomas F. El-Maraghi
% May 2004

if ~exist('sample')
   sample = 7.0/2.0;
end

% Determine the number of filter taps.
n = 2*round(sample * sigma)+1;

% Generate the x values.
x=1:n;
x=x-ceil(n/2);

% Sample the gaussian function to generate the filter taps.
g = exp(-(x.^2)/(2*sigma^2))/(sigma*sqrt(2*pi));

   