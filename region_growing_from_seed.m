
function j2 = region_growing_from_seed(im3, bin, reg_maxdist)
    % Perform region growing starting from seed points (centroids with area > 300)

    % Step 1: Binarize the image
    im4 = imbinarize(im3, bin);  % Binarize image with a threshold
    i22 = im2uint8(im4);  % Convert the binary image to uint8 for visualization
    
    % Step 2: Extract region properties (area and centroid)
    infobin = regionprops(im4, 'Area', 'Centroid');
    areabinary = cat(1, infobin.Area);  % Get the areas of the regions
    Centroidbinary = cat(1, infobin.Centroid);  % Get the centroids
    
    % Step 3: Filter regions with area greater than 300
    filtered_regions = Centroidbinary(areabinary > 300, :);
    final_grown_region = round(filtered_regions);

    % Step 4: Initialize output image
    j2 = zeros(size(im3));  % Create an empty matrix the same size as the input image
    
    % Step 5: Loop through filtered regions for region growing
    for lop = 1:size(final_grown_region, 1)
        x = final_grown_region(lop, 2);
        y = final_grown_region(lop, 1);
        J = region_grow(im3, x, y, reg_maxdist);
        j2 = j2 + J;  % Accumulate the grown regions
    end
end

function J = region_grow(I, x, y, reg_maxdist)
    % Grow region starting from (x, y) until intensity difference exceeds reg_maxdist

    if nargin < 4
        reg_maxdist = 0.2; % Default threshold
    end
    
    I = im2double(I);
    J = zeros(size(I)); % Output segmented region
    reg_mean = I(x, y); % The mean of the segmented region
    reg_size = 1; % Number of pixels in region
    neg_list = [x, y, I(x, y)]; % Neighbors list
    
    % Neighbor locations (4-connected)
    neigb = [-1 0; 1 0; 0 -1; 0 1];
    
    % Region growing loop
    while ~isempty(neg_list) && reg_size < numel(I)
        % Pop the first neighbor
        xn = neg_list(1, 1);
        yn = neg_list(1, 2);
        pix_value = neg_list(1, 3);
        neg_list(1, :) = []; % Remove the processed neighbor
        
        % Add the pixel to the region if not already added
        if J(xn, yn) == 0
            J(xn, yn) = 1;
            reg_size = reg_size + 1;
            reg_mean = (reg_mean * (reg_size - 1) + pix_value) / reg_size;
            
            % Add neighbors to the list
            for j = 1:4
                x_new = xn + neigb(j, 1);
                y_new = yn + neigb(j, 2);
                
                % Check if the neighbor is inside the image
                if x_new >= 1 && y_new >= 1 && x_new <= size(I, 1) && y_new <= size(I, 2)
                    if J(x_new, y_new) == 0 && abs(I(x_new, y_new) - reg_mean) < reg_maxdist
                        neg_list = [neg_list; x_new, y_new, I(x_new, y_new)];
                    end
                end
            end
        end
    end
end
