function [outputArg1, outputArg2, h, g, binLocations] = djusting_by_histogram_roi(inputArg1,pathg,pathh,pathbinLocations)
    % Function to process the input data and apply transformations to the image
    % This function takes an input image, processes it, and returns transformed data
    % Detailed explanation of the function is provided below.

    % Load external data from Excel files
%[pathg,pathh,pathbinLocations]=Extract_histogram_from_ROI();
    g = xlsread(pathg);
    h = xlsread(pathh);
    binLocations = xlsread(pathbinLocations);
    
    % Get the input image
    xw = inputArg1;
    
    % Get the dimensions of the image
    [ssii1, ssii2] = size(xw);
    
    % Flatten the image into a vector for processing
    ImgVector = reshape(xw, 1, []);
    ImgVector = im2uint8(ImgVector);  % Convert the image vector to uint8 type
    
    % Calculate the total number of elements in the image
    mult = ssii1 * ssii2;
    
    % Initialize a new image vector for transformed data
    ImgVector1 = uint8(zeros(1, mult));
    
    % Process each pixel in the image vector
    for i = 1:mult
        if ImgVector(1, i) < 156
            % If the pixel value is less than 156, apply the transformation based on `g`
            ind = ImgVector(1, i) + 1;
            ImgVector1(1, i) = fix((255 / g(108, 1)) * g(ind, 1));
        else
            % If the pixel value is greater or equal to 156, retain the original value
            ImgVector1(1, i) = ImgVector(1, i);
        end
    end
    
    % Convert the transformed vector back to an image matrix
    ImgVector2 = uint8(ImgVector1);
    
    % Normalize the image to the range [0, 1]
    ImgVector2 = mat2gray(reshape(ImgVector2, ssii1, []));
    
    % Set the output arguments
    outputArg1 = ImgVector2;
    outputArg2 = 1;  % This seems to be a placeholder; can be changed as needed
end
