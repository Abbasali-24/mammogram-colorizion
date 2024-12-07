
function [outputImage] = removePectoralMuscle1(inputImage)
    % Enhanced function for breast segmentation from a mammography image
    % It includes image preprocessing, breast segmentation, and mask generation.

%% 
    %% Image Preprocessing
    preprocessedImage = inputImage;
    if(size(preprocessedImage, 3) ~= 1)
        preprocessedImage = rgb2gray(preprocessedImage);
    end
    preprocessedImage = medfilt2(preprocessedImage, [3 3]); % Apply median filter with specified size

    % Convert to double and set parameters
    preprocessedImage = im2double(preprocessedImage); % Convert to double for precision
    isMLOView = true; % true for MLO view
    isFFDM = true; % true only for FFDM images

    %% Breast Segmentation
    % Segment boundary and chest wall
    [segmentationMask, ~, ~] = segmentBreastCombined(preprocessedImage, isMLOView, isFFDM);

    % Subtract pectoral muscle from original image
    outputImage = preprocessedImage .* segmentationMask;

    % Display the result
    figure;
    subplot(1, 2, 1);
    imshow(inputImage, []);
    title('Original Image');
    
    subplot(1, 2, 2);
    imshow(outputImage, []);
    title('Pectoral Muscle Removed');
end

function [segmentationMask, breastContour, chestWall] = segmentBreastCombined(inputMammogram, isMLOView, isFFDM)
    % Breast segmentation from mammography image, combines FFDM foreground
    % and other segmentation steps.
    
    % By default, assume that input is FFDM
    if nargin < 3
        isFFDM = true;
    end

    % Flip image if necessary
    isImageFlipped = isRightMammogram (inputMammogram);
    if isImageFlipped
        inputMammogram = fliplr(inputMammogram);
    end

    % Breast boundary detection
    if isFFDM
        [segmentationMask, breastContour] = detectBreastRegion (inputMammogram, isMLOView);
    else
        [segmentationMask, breastContour] = sfmForeground(inputMammogram, isMLOView);
    end

    breastContour.flip = isImageFlipped;

    % Breast chest wall detection
    if isMLOView
        [chestWallMask, chestWall] = segmentChestWall (inputMammogram, breastContour);
    else
        chestWallMask = true(size(segmentationMask));
        chestWall = struct('m', size(segmentationMask, 1), 'b', 0);
    end

    segmentationMask = segmentationMask & chestWallMask;

    % Flip back mask if necessary
    if isImageFlipped
        segmentationMask = fliplr(segmentationMask);
    end
end

function [segmentationMask, breastContour] = detectBreastRegion(inputMammogram, isMLOView)
    % Detect breast foreground in FFDM image with enhanced logic

    % Parameters:
    numBins = 1000; % Number of bins for image histogram

    if (nargin < 2) || isempty(isMLOView)
        isMLOView = false;
    end

    % Find intensity threshold
    warning off
    intensityThreshold = calculateIntensityThreshold (inputMammogram(:), numBins);
    warning on

    % Find mask
    initialMask = (inputMammogram >= max([intensityThreshold, 0]));

    % Remove artifacts and holes in the mask
    segmentationMask = refineSegmentationMask (initialMask);

    % Compute contour
    breastContour = extractBreastContour(segmentationMask, isMLOView);
    breastContour.th = intensityThreshold;
    segmentationMask(round(max(breastContour.y)):end, :) = false;
end

function intensityThreshold = calculateIntensityThreshold(x, numBins)
    % Calculate the lower threshold of the image intensity for segmentation

    % Minimum and maximum intensities:
    x_min = min(x);
    x_max = max(x);

    % Relative frequency:
    xi = linspace(x_min, x_max, numBins)';
    n = histc(x, xi);

    % Smooth histogram:
    n = conv(n, gausswin(25), 'same');
    n = n / max(n);

    % Find threshold by fitting a Gaussian to the histogram peak(s) below the mean
    xsup = min(mean(x), max(prctile(x, 30), .2));
    [~, ipeaks] = findpeaks(n .* (xi <= xsup), 'MinPeakHeight', .35);

    if (numel(ipeaks) == 1) % Only one peak
        select = (n > .35) & (xi < xsup);
        f = fit(xi(select), n(select), 'gauss1');
        intensityThreshold = f.b1 + sqrt(f.c1^2 * (log(f.a1) - log(0.05)));        

    elseif (numel(ipeaks) > 1) % Two peaks
        % Find minimum between peaks
        n_min = min(n(min(ipeaks):max(ipeaks)));
        i_min = find(n == n_min, 1);

        % Adjust second peak
        select = (xi >= xi(i_min)) & (n > .35) & (xi < xsup);
        f = fit(xi(select), n(select), 'gauss1');
        intensityThreshold = f.b1 + sqrt(f.c1^2 * (log(f.a1) - log(0.05)));

    elseif isempty(ipeaks) % No peaks
        n_max = max(n);
        i_th = find(n < 0.05 * n_max, 1);
        intensityThreshold = xi(i_th);
    end
end

function segmentationMask = refineSegmentationMask(initialMask)
    % Clean up the initial mask by removing artifacts and holes

    % Remove first and last row
    initialMask(1, :) = false;
    initialMask(end, :) = false;
    initialMask = imerode(initialMask, ones(5));

    % Keep biggest region
    cc = bwconncomp(initialMask);
    stats = regionprops(cc, 'area');
    idx = find([stats.Area] == max([stats.Area]), 1);
    initialMask = ismember(labelmatrix(cc), idx);

    % Remove spurious holes:
    initialMask = imdilate(initialMask, ones(5));
    initialMask = imclose(initialMask, ones(5));
    segmentationMask = imfill(initialMask, 'holes');
end

function flag = isRightMammogram(inputMammogram)
    % Determine whether mammogram is left or right.
    % TRUE if the mammogram is in the right side of the image.

    inputMammogram(inputMammogram > .95 * max(inputMammogram(:))) = 0;
    s = sum(inputMammogram);
    numSamples = round(0.5 * size(inputMammogram, 2));

    flag = sum(s(1:numSamples)) < sum(s(numSamples + 1:end));
end

function breastContour = extractBreastContour(segmentationMask, isMLOView)
    % Retrieve breast contour

    % Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K_th = 0.05; % Curvature threshold (higher is more strict)
    elim = 8;    % No of pixels to edge limit
    npts = 100;  % number of contour points to return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Find contour points:
    B = bwboundaries(segmentationMask);
    remov = (B{1}(:, 1) <= elim) | (B{1}(:, 2) <= elim) | (B{1}(:, 1) >= (size(segmentationMask, 1) - elim));

    ys = B{1}(~remov, 1);
    xs = B{1}(~remov, 2);

    % sub-sample and smooth contour:
    numSamples = length(xs);
    xs = smooth(interp1(xs, linspace(1, numSamples, npts), 'pchip'));
    ys = smooth(interp1(ys, linspace(1, numSamples, npts), 'pchip'));

    % Crop contour by curvature analysis
    [xc, yc, ycut] = cropBreastEdgeContour (xs, ys, K_th);

    if ~isMLOView
        breastContour.x = xs(:);
        breastContour.y = ys(:);
    else
        breastContour.x = xc(:);
        breastContour.y = yc(:);
    end

    breastContour.ycut = ycut;
    breastContour.size = size(segmentationMask);
end

function [xc, yc, ycut] = cropBreastEdgeContour(xs, ys, K_th, inputMammogram)
    % Crop breast contour using curvature thresholding

    dispflag = (nargin > 3);

    % Compute curvature k:
    dx1 = diff(xs);    % dx/dt
    dx2 = diff(dx1);   % d2x/dt2
    dy1 = diff(ys);    % dy/dt
    dy2 = diff(dy1);   % d2y/dt2
    dx1(1) = [];
    dy1(1) = [];

    k = (dx1 .* dy2 - dy1 .* dx2) ./ ((dx1 .^ 2 + dy1 .^ 2) .^ (1.5));

    % Cut contour points with curvature above the threshold K_th
    [kmin, i] = min(k);

    if (abs(kmin) > K_th) && (xs(i) < 0.4 * max(xs)) && (ys(i) > 0.5 * max(ys))
        ycut = floor(ys(i));
        xs(i + 1:end) = [];
        ys(i + 1:end) = [];
    else
        ycut = floor(max(ys));
    end    

    xc = xs;
    yc = ys;

    if dispflag
        figure
        subplot(121)
        imshow(mat2gray(inputMammogram)), hold on
        plot(xc, yc, 'g.')
        subplot(122)
        plot(k), title('Curvature'), hold on
        line([1 length(k)], [-K_th, -K_th], 'color', 'r')
    end
end

function [segmentationMask, chestWall] = segmentChestWall(inputMammogram, breastContour)

%Pre-process image:
processedImage = imdilate(inputMammogram, strel('disk',8,0));
processedImage = imfilter(processedImage, fspecial('gaussian',[5 5],1));

%%%% Detect pectoral line using Hough transform %%%%
%Find edges:
edge_map = edge(processedImage, 'canny', [], 2.0);



ymax = round(0.6*max(breastContour.y));
xmax = round(min(breastContour.x(breastContour.y<ymax)));
edge_map = edge_map(1:ymax, 1:xmax);

%remove lower diagonal:
[m,n] = size(edge_map);
[x, y] = meshgrid(1:n, 1:m);
yref = m - (m - 1)*(x-1)/(n-1);
edge_map(y>yref) = 0;


% Find edges coordinates:
[y, x] = find(edge_map);
%Quantize parameters space:
N = 128; %Number of quantization points
% rho_max = sqrt(sum(size(inputMammogram).^2));
rho_max = min(size(edge_map));
rho_min = 1;
theta_min = 20*pi/180;
theta_max = 45*pi/180;

%Compute accumulation array A
theta = linspace(theta_min, theta_max, N);
rho = linspace(rho_min, rho_max, N);
rho_k = x*cos(theta) + y*sin(theta);
A = histc(rho_k, rho);

%Find maximum in accumulation array:
[~,imax] = max(A(:));
[i,j] = ind2sub(size(A), imax);

%Get rect
T = theta(j);
R = (rho_max - rho_min)*(i - 1)/(N-1) + rho_min;
b = R/sin(T);
m = -cos(T)/sin(T);

%Get mask
[x, y] = meshgrid(1:size(inputMammogram,2), 1:size(inputMammogram,1));
segmentationMask = true(size(inputMammogram));
segmentationMask(y<b+m*x) = false;
chestWall.m = m;
chestWall.b = b;
end
