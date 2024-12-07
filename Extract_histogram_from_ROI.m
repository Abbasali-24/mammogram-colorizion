
function [file1, file2, file3] = Extract_histogram_from_ROI()

%% Extract histogram from ROI in mammogram dataset
% This code is used to extract the intensity histogram specifically for the tumor
% area in each mammogram that shows cancer in the dataset.


%% Read image data
h = double(zeros(256, 1));
g = double(zeros(256, 1));
 [path1] = get_excel_file_paths()
 
A = readtable(path1);

ihi = 0;
iihi = 0;

for i = 1:height(A)
    % Check if the image is marked as showing cancer
    if strncmpi(A.Status(i), 'Cancer', 1)
        % Skip images without a tumor contour
        if strncmpi(A.Tumour_Contour(i), '-', 1)
            ihi = ihi + 1;
        else
            % Read mammogram image and mask
            imagedir = char(fullfile('D:\', A.fullPath(i)));
            cancerimage = imread(imagedir);
            maskdir = char(fullfile('D:\', A.Tumour_Contour(i)));
            maskim = imread(maskdir);
            
            %% Adjust mask size to match image size
            [sz1, sz2] = size(cancerimage);
            [sz3, sz4] = size(maskim);
            zeroim = zeros(sz1, sz2);
            
            if sz1 * sz2 > sz3 * sz4
                zeroim(1:sz3, 1:sz4) = maskim;
            elseif sz1 * sz2 < sz3 * sz4
                zeroim = maskim(1:sz1, 1:sz2);
            else
                zeroim = maskim;
            end
            
            BW2 = imfill(zeroim, 'holes');
            
            %% Multiply mask with original image to extract tumor region
            Z = BW2 / 255;
            Z1 = immultiply(cancerimage, Z);
            
            % Calculate histogram of the tumor area
            [counts, binLocations] = imhist(Z1);
            counts(1) = 0; % Set background count to zero
            
            countsSum = sum(counts);
            countsNormalized = counts / countsSum;
            
            % Accumulate histograms across all tumor regions
            h = h + double(counts);
            g = g + countsNormalized;
        end
    end
end

% Adjust boundary values of histograms
h(1) = h(2);
h(256) = h(255);
g(1) = g(2);
g(256) = g(255);

% Plot histograms
% figure(1);
% bar(binLocations, h);
% title('Accumulated Histogram of Tumor Area (Raw Counts)');
% xlabel('Intensity Value');
% ylabel('Frequency');
% 
% figure(2);
% bar(binLocations, g);
% title('Accumulated Normalized Histogram of Tumor Area');
% xlabel('Intensity Value');
% ylabel('Normalized Frequency');

% Smooth histogram g
g1 = g * 0;
g1(1:122) = g(1:122);
for i = 122:255
    g1(i) = (g(i-1) + g(i) + g(i+1)) / 3;
end
file1='C:\Users\acer\Desktop\sh1.xlsx';
file2='C:\Users\acer\Desktop\sh2.xlsx';
file3='C:\Users\acer\Desktop\sh3.xlsx';

%Write histograms to Excel files 
xlswrite(file1, g);
xlswrite(file2, h);
xlswrite(file3, uint8(binLocations));

end
 