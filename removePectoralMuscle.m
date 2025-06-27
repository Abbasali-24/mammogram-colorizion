 
 function I_out = removePectoralMuscle(I)
% removePectoralMuscle - Removes the pectoral muscle from a mammogram image.
%
% INPUT:
%   I - Input grayscale mammogram image (uint8)
%
% OUTPUT:
%   I_out - Mammogram image with the pectoral muscle region removed (set to 0)

    % Convert to grayscale if needed
    if size(I, 3) == 3
        I = rgb2gray(I);
    end

    % Apply median filter to reduce noise
    I_filtered = medfilt2(I, [3 3]);

    % Apply global threshold to isolate bright regions (pectoral muscle)
    threshold_value = prctile(I_filtered(:), 95);
    BW = I_filtered > threshold_value;
imshow (BW)
    % Determine which side contains the pectoral muscle
    patch_size = 50;
    left_region_mean = mean2(I_filtered(1:patch_size, 1:patch_size));
    right_region_mean = mean2(I_filtered(1:patch_size, end-patch_size+1:end));

    if left_region_mean > right_region_mean
        seed_col = 10; % Muscle on left
    else
        seed_col = size(I, 2) - 10; % Muscle on right
    end
    seed_row = 10;

    % Use region growing to identify pectoral muscle
    muscle_region = bwselect(BW, seed_col, seed_row,8);

    % Remove pectoral muscle by zeroing out detected region
    I_out = I_filtered;
    I_out(muscle_region) = 0;
end
