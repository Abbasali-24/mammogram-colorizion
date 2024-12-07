function [outputArg1  ] = removecoulm (x0 )
 %% remove empty roe and colum

gray_image = rgb2gray(x0);
 row_sum = sum(gray_image, 2); % Sum along rows
column_sum = sum(gray_image, 1); % Sum along columns
non_empty_rows = find(row_sum > 1000);
non_empty_columns = find(column_sum > 1000);

cropped_image = x0(non_empty_rows, non_empty_columns, :);


outputArg1=cropped_image;
end

