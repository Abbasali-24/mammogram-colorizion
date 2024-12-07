
function [enhancedImage] = CONTRAST_ENHANCEMENT(inputImage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function applies an image enhancement algorithm using median-mean
% based sub-image clipped histogram equalization for contrast enhancement,
% inspired by the work of Kuldeep Singh, Rajiv Kapoor.
% enhancedImage: Enhanced output image with contrast improved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 1: Histogram Calculation
numGrayLevels = 256;  % Number of gray levels
[numRows, numCols] = size(inputImage);  % Dimensions of the input image
totalPixels = numRows * numCols;  % Total number of pixels
reshapedImage = reshape(inputImage, totalPixels, 1);  % Reshape image into a column vector
pdf = hist(reshapedImage, 0:numGrayLevels-1);  % Calculate the probability density function (PDF)

%% Step 2: Clipping Process (to limit extreme values)
clipThreshold = mean(pdf);  % Calculate the mean pixel count for gray levels
clipThreshold = round(clipThreshold);  % Round to nearest integer
clippedHist = zeros(1, numGrayLevels);  % Initialize intermediate histogram for clipping

for i = 1:numGrayLevels
    if pdf(i) > clipThreshold
        clippedHist(i) = clipThreshold;  % Clip histogram values above threshold
    else
        clippedHist(i) = pdf(i);  % Keep other values unchanged
    end
end

%% Step 3: Cumulative Distribution Function (CDF) for Histogram Equalization
cumulativeSum = zeros(numGrayLevels, 1);
totalpixels = sum(pdf, 2);
imageHist = imhist(inputImage);
cumulativeSum(1) = imageHist(1);

if (cumulativeSum(1) >= (round(totalpixels / 2)))
    medianValue = 0;
else
    for k = 2:numGrayLevels
        cumulativeSum(k) = cumulativeSum(k-1) + imageHist(k);
        if (cumulativeSum(k) >= round(totalpixels / 2))
            medianValue = k - 1;
            break;
        end
    end
end

%% Step 4: Sub-Histogram Equalization Using Median (medianValue)
[numRowsImage, numColsImage] = size(inputImage);
n_L1 = sum(imageHist(1:medianValue+1));
P_L1 = imageHist(1:medianValue+1) / n_L1;
n_U1 = sum(imageHist(medianValue+2:numGrayLevels));
P_U1 = imageHist(medianValue+2:numGrayLevels) / n_U1;
lowerCDF = zeros(medianValue+1, 1);  % Lower histogram CDF
upperCDF = zeros(numGrayLevels-medianValue-1, 1);  % Upper histogram CDF

% Apply sub-histogram equalization based on median (medianValue)
for u = 1:numRowsImage
    for t = 1:numColsImage 
        if (inputImage(u,t) < (medianValue+1))
            lowerCDF(1) = 0;
            for k = 2:(medianValue+1)
                lowerCDF(k) = lowerCDF(k-1) + (P_L1(k) * (k-1));
            end
        else
            upperCDF(1) = (medianValue+1) * P_U1(1);
            for k = 2:(numGrayLevels-medianValue-1)
                upperCDF(k) = upperCDF(k-1) + (P_U1(k) * (medianValue + k));
            end
        end
    end
end

lowerClip = round(lowerCDF(medianValue+1));
upperClip = round(upperCDF(numGrayLevels-medianValue-1));

%% Step 5: Final Equalization Based on Sub-Histograms
[imageWidth, imageHeight] = size(inputImage);
enhancedImage = zeros(size(inputImage));  

lowerClipCDF = zeros(1, lowerClip+1);
lowerClipHistSum = sum(clippedHist(1:lowerClip+1));
lowerClipPDF = clippedHist(1:lowerClip+1) / lowerClipHistSum;
lowerClipCDF(1) = lowerClipPDF(1);

for imageHist = 2:length(lowerClipPDF)
    lowerClipCDF(imageHist) = lowerClipPDF(imageHist) + lowerClipCDF(imageHist-1);
end

if medianValue ~= lowerClip
    upperClipCDF = zeros(1, (medianValue-lowerClip));
    upperClipHistSum = sum(clippedHist(lowerClip+2:medianValue));
    upperClipPDF = clippedHist(lowerClip+2:medianValue) / upperClipHistSum;
    upperClipCDF(1) = upperClipPDF(1);
    for imageHist = 2:(length(upperClipPDF))
        upperClipCDF(imageHist) = upperClipPDF(imageHist) + upperClipCDF(imageHist-1);
    end
end

lowerUnclipCDF = zeros(1, upperClip-medianValue);
upperUnclipCDF = zeros(1, (numGrayLevels-upperClip-1));
lowerUnclipHistSum = sum(clippedHist(medianValue+1:upperClip+1));
upperUnclipHistSum = sum(clippedHist(upperClip+2:numGrayLevels));
lowerUnclipPDF = clippedHist(medianValue+1:upperClip+1) / lowerUnclipHistSum;
upperUnclipPDF = clippedHist(upperClip+2:numGrayLevels) / upperUnclipHistSum;

lowerUnclipCDF(1) = lowerUnclipPDF(1);
for imageHist = 2:length(lowerUnclipPDF)
    lowerUnclipCDF(imageHist) = lowerUnclipPDF(imageHist) + lowerUnclipCDF(imageHist-1);
end

upperUnclipCDF(1) = upperUnclipPDF(1);
for imageHist = 2:(length(upperUnclipPDF))
    upperUnclipCDF(imageHist) = upperUnclipPDF(imageHist) + upperUnclipCDF(imageHist-1);
end

for imageHist = 1:imageWidth  % Apply the 4 equalized histograms 
    for s = 1:imageHeight
        if inputImage(imageHist,s) < (lowerClip+1)
            f = lowerClip * lowerClipCDF(inputImage(imageHist,s)+1);
            enhancedImage(imageHist,s) = round(f);
        elseif (inputImage(imageHist,s) >= (lowerClip+1)) && (inputImage(imageHist,s) < (medianValue+1))
            f = (lowerClip+1) + (medianValue-lowerClip+1) * upperClipCDF((inputImage(imageHist,s)-(lowerClip+1)) + 1);
            enhancedImage(imageHist,s) = round(f);
        elseif (inputImage(imageHist,s) >= medianValue+1) && (inputImage(imageHist,s) < (upperClip+1))
            f = (medianValue+1) + (upperClip-medianValue+1) * lowerUnclipCDF((inputImage(imageHist,s)-(medianValue+1)) + 1);
            enhancedImage(imageHist,s) = round(f);
        else
            f = (upperClip+1) + (numGrayLevels-upperClip+1) * upperUnclipCDF((inputImage(imageHist,s)-(upperClip+1)) + 1);
            enhancedImage(imageHist,s) = round(f);
        end
    end
end

%% Step 6: Output Type Adjustment (Ensure Output Matches Input Type)
switch class(inputImage)
    case 'uint8'
        enhancedImage = uint8(enhancedImage);
    case 'uint16'
        enhancedImage = uint16(enhancedImage);
    case 'int16'
        enhancedImage = int16(enhancedImage);
    case 'single'
        enhancedImage = single(enhancedImage);
    otherwise
        error('Unsupported input image class.');
end

end
