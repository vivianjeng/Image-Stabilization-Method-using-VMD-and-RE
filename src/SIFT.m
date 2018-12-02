function [ pos, scale, orient, desc ] = SIFT( im, octaves, intervals, object_mask, contrast_threshold, curvature_threshold, interactive, picNum )

% [ pos, scale, orient, desc ] = SIFT( im, octaves, intervals, object_mask, contrast_threshold, curvature_threshold, interactive )
%
% Apply David Lowe's Scale Invariant Feature Transform (SIFT) algorithm
% to a grayscale image.  This algorithm takes a grayscale image as input
% and returns a set of scale- and rotationally-invariant keypoints
% allong with their corresponding feature descriptors.
%
% This implementation is based on:
% [1] David G. Lowe, "Distinctive Image Features from Sacle-Invariant Keypoints",
%     accepted for publicatoin in the International Journal of Computer
%     Vision, 2004.
% [2] David G. Lowe, "Object Recognition from Local Scale-Invariant Features",
%     Proc. of the International Conference on Computer Vision, Corfu,
%     September 1999.
%
% Input:
% im - the input image, with pixel values normalize to lie betwen [0,1].
% octaves - the number of octaves to search for keypoints (default = 4).  
% intervals - the number of geometrically sampled intervals to divide
%   each octave into when searching for keypoints (default = 2).
% object_mask - a binary mask specifying the location of the object in
%   the image to search for keypoints on.  If not specified, the whole
%   image is searched.
% contrast_threshold - the threshold on the contrast of the DOG extrema
%   before classifying them as keypoints (default = 0.03).
% curvature_threshold - the upper bound on the ratio between the principal
%   curvatures of the DOG extrema before classifying it as a keypoint
%   (default = 10.0).
% interactive - >= 1 displays progress and timing information,
%   >= 2 displays intermediate results of the algorihtm (default = 1).
%
% Output:
% pos - an Nx2 matrix containing the (x,y) coordinates of the keypoints
%   stored in rows.
% scale - an Nx3 matrix with rows describing the scale of each keypoint (i.e.,
%   first column specifies the octave, second column specifies the interval, and
%   third column specifies sigma).
% orient - a Nx1 vector containing the orientations of the keypoints [-pi,pi).
% desc - an Nx128 matrix with rows containing the feature descriptors 
%   corresponding to the keypoints.
%
% Thomas F. El-Maraghi
% May 2004

% assign default values to the input variables
if ~exist('octaves')
   octaves = 4;
end
if ~exist('intervals')
   intervals = 2;
end
if ~exist('object_mask')
   object_mask = ones(size(im));
end
if size(object_mask) ~= size(im)
   object_mask = ones(size(im));
end
if ~exist('contrast_threshold')
   contrast_threshold = 0.02;
end
if ~exist('curvature_threshold')
   curvature_threshold = 10.0;
end
if ~exist('interactive')
   interactive = 1;
end

% Check that the image is normalized to [0,1]
if( (min(im(:)) < 0) | (max(im(:)) > 1) )
   fprintf( 2, 'Warning: image not normalized to [0,1].\n' );
end

% Blur the image with a standard deviation of 0.5 to prevent aliasing
% and then upsample the image by a factor of 2 using linear interpolation.
% Lowe claims that this increases the number of stable keypoints by 
% a factor of 4.
if interactive >= 1
   fprintf( 2, 'Doubling image size for first octave...\n' );
end
tic;
antialias_sigma = 0.5;
if antialias_sigma == 0
   signal = im;
else
   g = gaussian_filter( antialias_sigma );
   if exist('corrsep') == 3
	   signal = corrsep( g, g, im );
   else
      signal = conv2( g, g, im, 'same' );
   end
end
signal = im;
[X Y] = meshgrid( 1:0.5:size(signal,2), 1:0.5:size(signal,1) );
signal = interp2( signal, X, Y, '*linear' );   
subsample = [0.5]; % subsampling rate for doubled image is 1/2


% The next step of the algorithm is to generate the gaussian and difference-of-
% gaussian (DOG) pyramids.  These pyramids will be stored as two cell arrays,
% gauss_pyr{orient,interval} and DOG_pyr{orient,interval}, respectively.  In order
% to detect keypoints on s intervals per octave, we must generate s+3 blurred
% images in the gaussian pyramid.  This is becuase s+3 blurred images generates
% s+2 DOG images, and two images are needed (one at the highest and one lowest scales 
% of the octave) for extrema detection.

% Generate the first image of the first octave of the gaussian pyramid
% by preblurring the doubled image with a gaussian with a standard deviation
% of 1.6.  This choice for sigma is a trade off between repeatability and
% efficiency.
if interactive >= 1
   fprintf( 2, 'Prebluring image...\n' );
end
preblur_sigma = sqrt(sqrt(2)^2 - (2*antialias_sigma)^2);
if preblur_sigma == 0
   gauss_pyr{1,1} = signal;
else
   g = gaussian_filter( preblur_sigma );
   if exist('corrsep') == 3
      gauss_pyr{1,1} = corrsep( g, g, signal );
   else
      gauss_pyr{1,1} = conv2( g, g, signal, 'same' );
   end
end
clear signal
pre_time = toc;
if interactive >= 1
   fprintf( 2, 'Preprocessing time %.2f seconds.\n', pre_time );
end

% The initial blurring for the first image of the first octave of the pyramid.
initial_sigma = sqrt( (2*antialias_sigma)^2 + preblur_sigma^2 );

% Keep track of the absolute sigma for the octave and scale
absolute_sigma = zeros(octaves,intervals+3);
absolute_sigma(1,1) = initial_sigma * subsample(1);

% Keep track of the filter sizes and standard deviations used to generate the pyramid
filter_size = zeros(octaves,intervals+3);
filter_sigma = zeros(octaves,intervals+3);

% Generate the remaining levels of the geometrically sampled gaussian and DOG pyramids
if interactive >= 1
   fprintf( 2, 'Expanding the Gaussian and DOG pyramids...\n' );
end
tic;
for octave = 1:octaves
   if interactive >= 1
      fprintf( 2, '\tProcessing octave %d: image size %d x %d subsample %.1f\n', octave, size(gauss_pyr{octave,1},2), size(gauss_pyr{octave,1},1), subsample(octave) );
      fprintf( 2, '\t\tInterval 1 sigma %f\n', absolute_sigma(octave,1) );
   end   
   sigma = initial_sigma;
   g = gaussian_filter( sigma );
   filter_size( octave, 1 ) = length(g);
   filter_sigma( octave, 1 ) = sigma;
   DOG_pyr{octave} = zeros(size(gauss_pyr{octave,1},1),size(gauss_pyr{octave,1},2),intervals+2);
   for interval = 2:(intervals+3)
      
      % Compute the standard deviation of the gaussian filter needed to produce the 
      % next level of the geometrically sampled pyramid.  Here, sigma_i+1 = k*sigma.
      % By definition of successive convolution, the required blurring sigma_f to
      % produce sigma_i+1 from sigma_i is:
      %
      %    sigma_i+1^2 = sigma_f,i^2 + sigma_i^2
      %  (k*sigma_i)^2 = sigma_f,i^2 + sigma_i^2
      %
      % therefore:
      %
      %      sigma_f,i = sqrt(k^2 - 1)sigma_i
      % 
      % where k = 2^(1/intervals) to span the octave, so:
      %
      %  sigma_f,i = sqrt(2^(2/intervals) - 1)sigma_i
      sigma_f = sqrt(2^(2/intervals) - 1)*sigma;
      g = gaussian_filter( sigma_f );
      sigma = (2^(1/intervals))*sigma;
      
      % Keep track of the absolute sigma
      absolute_sigma(octave,interval) = sigma * subsample(octave);
      
      % Store the size and standard deviation of the filter for later use
      filter_size(octave,interval) = length(g);
      filter_sigma(octave,interval) = sigma;
      
      if exist('corrsep') == 3
         gauss_pyr{octave,interval} = corrsep( g, g, gauss_pyr{octave,interval-1} );
      else
         gauss_pyr{octave,interval} = conv2( g, g, gauss_pyr{octave,interval-1}, 'same' );
      end      
      DOG_pyr{octave}(:,:,interval-1) = gauss_pyr{octave,interval} - gauss_pyr{octave,interval-1};
      
      if interactive >= 1
         fprintf( 2, '\t\tInterval %d sigma %f\n', interval, absolute_sigma(octave,interval) );
      end              
   end      
   if octave < octaves
      % The gaussian image 2 images from the top of the stack for
      % this octave have be blurred by 2*sigma.  Subsample this image by a 
      % factor of 2 to procuduce the first image of the next octave.
      sz = size(gauss_pyr{octave,intervals+1});
      [X Y] = meshgrid( 1:2:sz(2), 1:2:sz(1) );
      gauss_pyr{octave+1,1} = interp2(gauss_pyr{octave,intervals+1},X,Y,'*nearest'); 
      absolute_sigma(octave+1,1) = absolute_sigma(octave,intervals+1);
      subsample = [subsample subsample(end)*2];
   end      
end
pyr_time = toc;
if interactive >= 1
   fprintf( 2, 'Pryamid processing time %.2f seconds.\n', pyr_time );
end

% Display the gaussian pyramid when in interactive mode
if interactive >= 2
   sz = zeros(1,2);
   sz(2) = (intervals+3)*size(gauss_pyr{1,1},2);
   for octave = 1:octaves
      sz(1) = sz(1) + size(gauss_pyr{octave,1},1);
   end
   pic = zeros(sz);
   y = 1;
   for octave = 1:octaves
      x = 1;
      sz = size(gauss_pyr{octave,1});
      for interval = 1:(intervals + 3)
			pic(y:(y+sz(1)-1),x:(x+sz(2)-1)) = gauss_pyr{octave,interval};		         
         x = x + sz(2);
      end
      y = y + sz(1);
   end
   fig = figure;
   clf;
   showIm(pic);
   resizeImageFig( fig, size(pic), 0.25 );
   fprintf( 2, 'The gaussian pyramid (0.25 scale).\nPress any key to continue.\n' );
   %pause;
   close(fig)
end

% Display the DOG pyramid when in interactive mode
if interactive >= 2
   sz = zeros(1,2);
   sz(2) = (intervals+2)*size(DOG_pyr{1}(:,:,1),2);
   for octave = 1:octaves
      sz(1) = sz(1) + size(DOG_pyr{octave}(:,:,1),1);
   end
   pic = zeros(sz);
   y = 1;
   for octave = 1:octaves
      x = 1;
      sz = size(DOG_pyr{octave}(:,:,1));
      for interval = 1:(intervals + 2)
			pic(y:(y+sz(1)-1),x:(x+sz(2)-1)) = DOG_pyr{octave}(:,:,interval);		         
         x = x + sz(2);
      end
      y = y + sz(1);
   end
   fig = figure;
   clf;
   showIm(pic);
   resizeImageFig( fig, size(pic), 0.25 );
   fprintf( 2, 'The DOG pyramid (0.25 scale).\nPress any key to continue.\n' );
   %pause;
   close(fig)
end

% The next step is to detect local maxima in the DOG pyramid.  When
% a maximum is found, two tests are applied before labeling it as a 
% keypoint.  First, it must have sufficient contrast.  Second, it should
% not be and edge point (i.e., the ratio of principal curvatures at the
% extremum should be below a threshold).

% Compute threshold for the ratio of principle curvature test applied to
% the DOG extrema before classifying them as keypoints.
curvature_threshold = ((curvature_threshold + 1)^2)/curvature_threshold;

% 2nd derivative kernels 
xx = [ 1 -2  1 ];
yy = xx';
xy = [ 1 0 -1; 0 0 0; -1 0 1 ]/4;

% Coordinates of keypoints after each stage of processing for display
% in interactive mode.
raw_keypoints = [];
contrast_keypoints = [];
curve_keypoints = [];

% Detect local maxima in the DOG pyramid
if interactive >= 1
   fprintf( 2, 'Locating keypoints...\n' );
end
tic;
loc = cell(size(DOG_pyr)); % boolean maps of keypoints
for octave = 1:octaves
   if interactive >= 1
      fprintf( 2, '\tProcessing octave %d\n', octave );
   end
   for interval = 2:(intervals+1)
      keypoint_count = 0;
      contrast_mask = abs(DOG_pyr{octave}(:,:,interval)) >= contrast_threshold;
      loc{octave,interval} = zeros(size(DOG_pyr{octave}(:,:,interval)));
      if exist('corrsep') == 3
         edge = 1;
      else         
         edge = ceil(filter_size(octave,interval)/2);
      end      
      for y=(1+edge):(size(DOG_pyr{octave}(:,:,interval),1)-edge)         
         for x=(1+edge):(size(DOG_pyr{octave}(:,:,interval),2)-edge)
            % Only check for extrema where the object mask is 1
            if object_mask(round(y*subsample(octave)),round(x*subsample(octave))) == 1 
               
               % When not displaying intermediate results, perform the check that the current location
               % in the DOG pyramid is above the contrast threshold before checking
               % for an extrema for efficiency reasons.  Note: we could not make this
               % change of order if we were interpolating the locations of the extrema.
               if( (interactive >= 2) | (contrast_mask(y,x) == 1) ) 
                  
                  % Check for a max or a min across space and scale
                  tmp = DOG_pyr{octave}((y-1):(y+1),(x-1):(x+1),(interval-1):(interval+1));  
                  pt_val = tmp(2,2,2);
                  if( (pt_val == min(tmp(:))) | (pt_val == max(tmp(:))) )
                     % The point is a local extrema of the DOG image.  Store its coordinates for
                     % displaying keypoint location in interactive mode.
                     raw_keypoints = [raw_keypoints; x*subsample(octave) y*subsample(octave)];
                     if abs(DOG_pyr{octave}(y,x,interval)) >= contrast_threshold
                        % The DOG image at the extrema is above the contrast threshold.  Store 
                        % its coordinates for displaying keypoint locations in interactive mode.
                        contrast_keypoints = [contrast_keypoints; raw_keypoints(end,:)];
                        % Compute the entries of the Hessian matrix at the extrema location.
                        Dxx = sum(DOG_pyr{octave}(y,x-1:x+1,interval) .* xx);
                        Dyy = sum(DOG_pyr{octave}(y-1:y+1,x,interval) .* yy);
                        Dxy = sum(sum(DOG_pyr{octave}(y-1:y+1,x-1:x+1,interval) .* xy));
                        
                        % Compute the trace and the determinant of the Hessian.
                        Tr_H = Dxx + Dyy;
                        Det_H = Dxx*Dyy - Dxy^2;
                        
                        % Compute the ratio of the principal curvatures.
                        curvature_ratio = (Tr_H^2)/Det_H;
                        
                        if ((Det_H >= 0) & (curvature_ratio < curvature_threshold))
                           % The ratio of principal curvatures is below the threshold (i.e.,
                           % it is not an edge point).  Store its coordianates for displaying
                           % keypoint locations in interactive mode.
                           curve_keypoints = [curve_keypoints; raw_keypoints(end,:)];
                           % Set the loc map to 1 to at this point to indicate a keypoint.
                           loc{octave,interval}(y,x) = 1;
                           keypoint_count = keypoint_count + 1;
                        end
                     end                  
                  end
               end               
            end
         end         
      end
      if interactive >= 1
         fprintf( 2, '\t\t%d keypoints found on interval %d\n', keypoint_count, interval );
      end
   end
end
keypoint_time = toc;
if interactive >= 1
   fprintf( 2, 'Keypoint location time %.2f seconds.\n', keypoint_time );
end   

% Display results of extrema detection and keypoint filtering in interactive mode.
if interactive >= 2
   fig = figure();
   clf;
   showIm(im);
   hold on;
   plot(raw_keypoints(:,1),raw_keypoints(:,2),'y+');
   %saveas(fig, ['out_' int2str(picNum) '.jpg' ]);
   resizeImageFig( fig, size(im), 2 );
   fprintf( 2, 'DOG extrema (2x scale).\nPress any key to continue.\n' );
   %pause;
   close(fig);
   fig = figure;
   clf;
   showIm(im);
   hold on;
   plot(contrast_keypoints(:,1),contrast_keypoints(:,2),'y+');
   resizeImageFig( fig, size(im), 2 );
   fprintf( 2, 'Keypoints after removing low contrast extrema (2x scale).\nPress any key to continue.\n' );
   %pause;
   close(fig);
   fig = figure;
   clf;
   showIm(im);
   hold on;
   plot(curve_keypoints(:,1),curve_keypoints(:,2),'y+');
   resizeImageFig( fig, size(im), 2 );
   fprintf( 2, 'Keypoints after removing edge points using principal curvature filtering (2x scale).\nPress any key to continue.\n' );
   %pause;
   close(fig);  
end
clear raw_keypoints contrast_keypoints curve_keypoints

% The next step of the algorithm is to assign orientations to the keypoints.  For this,
% we histogram the gradient orientation over a region about each keypoint.
g = gaussian_filter( 1.5 * absolute_sigma(1,intervals+3) / subsample(1) );
zero_pad = ceil( length(g) / 2 );

% Compute the gradient direction and magnitude of the gaussian pyramid images
if interactive >= 1
   fprintf( 2, 'Computing gradient magnitude and orientation...\n' );
end
tic;
mag_thresh = zeros(size(gauss_pyr));
mag_pyr = cell(size(gauss_pyr));
grad_pyr = cell(size(gauss_pyr));
for octave = 1:octaves
   for interval = 2:(intervals+1)      
      % Compute x and y derivatives using pixel differences
      diff_x = 0.5*(gauss_pyr{octave,interval}(2:(end-1),3:(end))-gauss_pyr{octave,interval}(2:(end-1),1:(end-2)));
      diff_y = 0.5*(gauss_pyr{octave,interval}(3:(end),2:(end-1))-gauss_pyr{octave,interval}(1:(end-2),2:(end-1)));
      
      % Compute the magnitude of the gradient
      mag = zeros(size(gauss_pyr{octave,interval}));      
      mag(2:(end-1),2:(end-1)) = sqrt( diff_x .^ 2 + diff_y .^ 2 );
      
      % Store the magnitude of the gradient in the pyramid with zero padding
      mag_pyr{octave,interval} = zeros(size(mag)+2*zero_pad);
      mag_pyr{octave,interval}((zero_pad+1):(end-zero_pad),(zero_pad+1):(end-zero_pad)) = mag;      
      
      % Compute the orientation of the gradient
      grad = zeros(size(gauss_pyr{octave,interval}));
      grad(2:(end-1),2:(end-1)) = atan2( diff_y, diff_x );
      grad(find(grad == pi)) = -pi;
      
      % Store the orientation of the gradient in the pyramid with zero padding
      grad_pyr{octave,interval} = zeros(size(grad)+2*zero_pad);
      grad_pyr{octave,interval}((zero_pad+1):(end-zero_pad),(zero_pad+1):(end-zero_pad)) = grad;
   end
end
clear mag grad
grad_time = toc;
if interactive >= 1
   fprintf( 2, 'Gradient calculation time %.2f seconds.\n', grad_time );
end

% The next step of the algorithm is to assign orientations to the keypoints
% that have been located.  This is done by looking for peaks in histograms of
% gradient orientations in regions surrounding each keypoint.  A keypoint may be 
% assigned more than one orientation.  If it is, then two identical descriptors 
% are added to the database with different orientations.

% Set up the histogram bin centers for a 36 bin histogram.
num_bins = 36;
hist_step = 2*pi/num_bins;
hist_orient = [-pi:hist_step:(pi-hist_step)];

% Initialize the positions, orientations, and scale information
% of the keypoints to emtpy matrices.
pos = [];
orient = [];
scale = [];

% Assign orientations to the keypoints.
if interactive >= 1
   fprintf( 2, 'Assigining keypoint orientations...\n' );
end
tic;
for octave = 1:octaves
   if interactive >= 1
      fprintf( 2, '\tProcessing octave %d\n', octave );
   end
   for interval = 2:(intervals + 1)
      if interactive >= 1
         fprintf( 2, '\t\tProcessing interval %d ', interval );
      end            
      keypoint_count = 0;
      
      % Create a gaussian weighting mask with a standard deviation of 1/2 of
      % the filter size used to generate this level of the pyramid.
      g = gaussian_filter( 1.5 * absolute_sigma(octave,interval)/subsample(octave) );
      hf_sz = floor(length(g)/2);
      g = g'*g;      
      
      % Zero pad the keypoint location map.
      loc_pad = zeros(size(loc{octave,interval})+2*zero_pad);
      loc_pad((zero_pad+1):(end-zero_pad),(zero_pad+1):(end-zero_pad)) = loc{octave,interval};
      
      % Iterate over all the keypoints at this octave and orientation.
      [iy ix]=find(loc_pad==1);
      for k = 1:length(iy)
         % Histogram the gradient orientations for this keypoint weighted by the
         % gradient magnitude and the gaussian weighting mask.
         x = ix(k);
         y = iy(k);
         wght = g.*mag_pyr{octave,interval}((y-hf_sz):(y+hf_sz),(x-hf_sz):(x+hf_sz));
         grad_window = grad_pyr{octave,interval}((y-hf_sz):(y+hf_sz),(x-hf_sz):(x+hf_sz));
         orient_hist=zeros(length(hist_orient),1);
         for bin=1:length(hist_orient)
            % Compute the diference of the orientations mod pi
            diff = mod( grad_window - hist_orient(bin) + pi, 2*pi ) - pi;
            
            % Accumulate the histogram bins
            orient_hist(bin)=orient_hist(bin)+sum(sum(wght.*max(1 - abs(diff)/hist_step,0)));
            orient_hist(bin)=orient_hist(bin)+sum(sum(wght.*(abs(diff) <= hist_step)));
         end
         
         % Find peaks in the orientation histogram using nonmax suppression.
         peaks = orient_hist;        
         rot_right = [ peaks(end); peaks(1:end-1) ];
         rot_left = [ peaks(2:end); peaks(1) ];         
         peaks( find(peaks < rot_right) ) = 0;
         peaks( find(peaks < rot_left) ) = 0;
         
         % Extract the value and index of the largest peak. 
         [max_peak_val ipeak] = max(peaks);
         
         % Iterate over all peaks within 80% of the largest peak and add keypoints with
         % the orientation corresponding to those peaks to the keypoint list.
         peak_val = max_peak_val;
         while( peak_val > 0.8*max_peak_val )
            % Interpolate the peak by fitting a parabola to the three histogram values
            % closest to each peak.				            
            A = [];
            b = [];
            for j = -1:1
               A = [A; (hist_orient(ipeak)+hist_step*j).^2 (hist_orient(ipeak)+hist_step*j) 1];
	            bin = mod( ipeak + j + num_bins - 1, num_bins ) + 1;
               b = [b; orient_hist(bin)];
            end
            c = pinv(A)*b;
            max_orient = -c(2)/(2*c(1));
            while( max_orient < -pi )
               max_orient = max_orient + 2*pi;
            end
            while( max_orient >= pi )
               max_orient = max_orient - 2*pi;
            end            
            
            % Store the keypoint position, orientation, and scale information
            pos = [pos; [(x-zero_pad) (y-zero_pad)]*subsample(octave) ];
            orient = [orient; max_orient];
            scale = [scale; octave interval absolute_sigma(octave,interval)];
            keypoint_count = keypoint_count + 1;
            
            % Get the next peak
            peaks(ipeak) = 0;
            [peak_val ipeak] = max(peaks);
         end         
      end
      if interactive >= 1
         fprintf( 2, '(%d keypoints)\n', keypoint_count );
      end            
   end
end
clear loc loc_pad 
orient_time = toc;
if interactive >= 1
   fprintf( 2, 'Orientation assignment time %.2f seconds.\n', orient_time );
end

% Display the keypoints with scale and orientation in interactive mode.
if interactive >= 2
   fig = figure;
   clf;
   showIm(im);
   hold on;
   display_keypoints( pos, scale(:,3), orient, 'y' );
   resizeImageFig( fig, size(im), 2 );
   fprintf( 2, 'Final keypoints with scale and orientation (2x scale).\nPress any key to continue.\n' );
   %pause;
   close(fig);
end

% The final of the SIFT algorithm is to extract feature descriptors for the keypoints.
% The descriptors are a grid of gradient orientation histograms, where the sampling
% grid for the histograms is rotated to the main orientation of each keypoint.  The
% grid is a 4x4 array of 4x4 sample cells of 8 bin orientation histograms.  This 
% procduces 128 dimensional feature vectors.

% The orientation histograms have 8 bins
orient_bin_spacing = pi/4;
orient_angles = [-pi:orient_bin_spacing:(pi-orient_bin_spacing)];

% The feature grid is has 4x4 cells - feat_grid describes the cell center positions
grid_spacing = 4;
[x_coords y_coords] = meshgrid( [-6:grid_spacing:6] );
feat_grid = [x_coords(:) y_coords(:)]';
[x_coords y_coords] = meshgrid( [-(2*grid_spacing-0.5):(2*grid_spacing-0.5)] );
feat_samples = [x_coords(:) y_coords(:)]';
feat_window = 2*grid_spacing;

% Initialize the descriptor list to the empty matrix.
desc = [];

% Loop over all of the keypoints.
if interactive >= 1
   fprintf( 2, 'Computing keypoint feature descriptors for %d keypoints', size(pos,1) );
end
for k = 1:size(pos,1)
   x = pos(k,1)/subsample(scale(k,1));
   y = pos(k,2)/subsample(scale(k,1));   
   
   % Rotate the grid coordinates.
   M = [cos(orient(k)) -sin(orient(k)); sin(orient(k)) cos(orient(k))];
   feat_rot_grid = M*feat_grid + repmat([x; y],1,size(feat_grid,2));
   feat_rot_samples = M*feat_samples + repmat([x; y],1,size(feat_samples,2));
   
   % Initialize the feature descriptor.
   feat_desc = zeros(1,128);
   
   % Histogram the gradient orientation samples weighted by the gradient magnitude and
   % a gaussian with a standard deviation of 1/2 the feature window.  To avoid boundary
   % effects, each sample is accumulated into neighbouring bins weighted by 1-d in
   % all dimensions, where d is the distance from the center of the bin measured in
   % units of bin spacing.
   for s = 1:size(feat_rot_samples,2)
      x_sample = feat_rot_samples(1,s);
      y_sample = feat_rot_samples(2,s);
      
      % Interpolate the gradient at the sample position
      [X Y] = meshgrid( (x_sample-1):(x_sample+1), (y_sample-1):(y_sample+1) );
      G = interp2( gauss_pyr{scale(k,1),scale(k,2)}, X, Y, '*linear' );
      G(find(isnan(G))) = 0;
      diff_x = 0.5*(G(2,3) - G(2,1));
      diff_y = 0.5*(G(3,2) - G(1,2));
      mag_sample = sqrt( diff_x^2 + diff_y^2 );
      grad_sample = atan2( diff_y, diff_x );
      if grad_sample == pi
         grad_sample = -pi;
      end      
      
      % Compute the weighting for the x and y dimensions.
      x_wght = max(1 - (abs(feat_rot_grid(1,:) - x_sample)/grid_spacing), 0);
      y_wght = max(1 - (abs(feat_rot_grid(2,:) - y_sample)/grid_spacing), 0); 
      pos_wght = reshape(repmat(x_wght.*y_wght,8,1),1,128);
      
      % Compute the weighting for the orientation, rotating the gradient to the
      % main orientation to of the keypoint first, and then computing the difference
      % in angle to the histogram bin mod pi.
      diff = mod( grad_sample - orient(k) - orient_angles + pi, 2*pi ) - pi;
      orient_wght = max(1 - abs(diff)/orient_bin_spacing,0);
      orient_wght = repmat(orient_wght,1,16);         
      
      % Compute the gaussian weighting.
      g = exp(-((x_sample-x)^2+(y_sample-y)^2)/(2*feat_window^2))/(2*pi*feat_window^2);
      
      % Accumulate the histogram bins.
      feat_desc = feat_desc + pos_wght.*orient_wght*g*mag_sample;
   end
   
   % Normalize the feature descriptor to a unit vector to make the descriptor invariant
   % to affine changes in illumination.
   feat_desc = feat_desc / norm(feat_desc);
   
   % Threshold the large components in the descriptor to 0.2 and then renormalize
   % to reduce the influence of large gradient magnitudes on the descriptor.
   feat_desc( find(feat_desc > 0.2) ) = 0.2;
   feat_desc = feat_desc / norm(feat_desc);
   
   % Store the descriptor.
   desc = [desc; feat_desc];
   if (interactive >= 1) & (mod(k,25) == 0)
      fprintf( 2, '.' );
   end
end
desc_time = toc;

% Adjust for the sample offset
sample_offset = -(subsample - 1);
for k = 1:size(pos,1)
   pos(k,:) = pos(k,:) + sample_offset(scale(k,1));
end

% Return only the absolute scale
if size(pos,1) > 0
	scale = scale(:,3);
end
   
% Display summary in interactive mode.
if interactive >= 1
   fprintf( 2, '\nDescriptor processing time %.2f seconds.\n', desc_time );
   fprintf( 2, 'Processing time summary:\n' );
   fprintf( 2, '\tPreprocessing:\t%.2f s\n', pre_time );
   fprintf( 2, '\tPyramid:\t%.2f s\n', pyr_time );
   fprintf( 2, '\tKeypoints:\t%.2f s\n', keypoint_time );
   fprintf( 2, '\tGradient:\t%.2f s\n', grad_time );
   fprintf( 2, '\tOrientation:\t%.2f s\n', orient_time );
   fprintf( 2, '\tDescriptor:\t%.2f s\n', desc_time );
   fprintf( 2, 'Total processing time %.2f seconds.\n', pre_time + pyr_time + keypoint_time + grad_time + orient_time + desc_time );
end

