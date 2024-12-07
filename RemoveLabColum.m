function [outputArg1] = RemoveLabColum(x0)
% REMOVE LAB AND EMPTY ROWS/COLUMNS FROM IMAGE
%   Removes the label and any empty rows or columns from the input image.

% Remove the label
x01 = imbinarize(x0, 0.1);
info = regionprops(x01);
a = cat(1, info.Area);
[~, linfo] = max(a);
x011 = info(linfo).Centroid;
x2 = bwselect(x01, x011(1), x011(2), 8);
x0 = immultiply(x0, x2);
x0(1, :) = 0;

% Remove empty rows and columns
[siz1, siz2] = size(x0);
tst1 = zeros(siz1, siz2);

% Detect empty rows
r1 = find(all(x0 == tst1, 2));

% Detect empty columns
r2 = find(all(x0 == tst1, 1));

% Delete empty rows and columns
x0(:, r2) = [];
x0(r1, :) = [];

% Output the processed image
outputArg1 = x0;
end
