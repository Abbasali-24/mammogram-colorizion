clc;  % Clear the command window.
close all;  % Close all figures (except those of imtool).
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;
%     setup;
%% get_excel_file_paths
% path 1 this is the execel dataset filel  for image
% path 2  to write output color  image

% Initialize output variables
[path1]=get_excel_file_paths; %this is the execel path  dataset filel  for image
[~,output_folder] = get_excel_file_paths; %to write output color  image

%this is the path to path of histogram that calclated from roi used in djusting_by_histogram_roi function
[pathg, pathh, pathbinLocations] = Extract_histogram_from_ROI()

A=readtable(path1);   % this is the execel dataset filel  for image
lop = [];
lo = 1;
for i1 = 1:height(A) % this loop from 1 : end dataset fiel
    %% defintion the brast image direction and it information in data set like mask image ..
    % Define image directory from dataset
    imagedir = char(fullfile('D:\', A.fullPath(i1)));
    breastrimage = imread(imagedir);
    [sa1, sa2] = size(imagedir);
    % Define mask  directory from dataset for each mamogram image
    maskfromimag = char(imagedir(1:(sa2-4)));
    typeimage=char(imagedir(4));
    maskfromimag1 = ['_Mask.jpg'];
    maskfromimag2 = [maskfromimag, maskfromimag1];
    
    if exist(maskfromimag2, 'file')
        % Load mask image
        maskim = imread(maskfromimag2);
        
        
        %% Process the cancer image
        %% remove lab and empty row and coulum removelabcoulum and noise  and Pectoral Muscle1
        
        imagenumer = breastrimage;
        [sizce1, sizcen2] = size(breastrimage);
        imagenumer(1:100, :) = [];
        imagenumer((sizce1-100):end, :) = [];
        
        %% Apply various image processing steps
        imagenumer = medfilt2(imagenumer, [3 3]) 
        im1 = RemoveLabColum(imagenumer); %% remove lab and empty row and coloum
        im2= removePectoralMuscle1(im1);%% remove pectoral muscles
        im2 = im2uint8(im2);
        %% Enhance the perceptibility of mammography images
       contrasrenhance=CONTRAST_ENHANCEMENT(im2);

        K = medfilt2(contrasrenhance, [8 8]);
        I = im2double(K);
        I = wiener2(I, [5 5]);
        qqq = (K - contrasrenhance);
        im3 = im2uint8(qqq + contrasrenhance);

window_size = round(size(im3, 2) / 8);  % 1/8 of image width
T = adaptthresh(I, 0.15, 'NeighborhoodSize', window_size, 'Statistic', 'mean');
 
imshow(contrasrenhance)
        %% Perform region growing and combine images
        %Select suspicious lesions in mammograms by thresholding and the flood-fill technique
        [final_image] = region_growing_from_seed(im3,T,0.2);
        
        im4=im2uint8(final_image*0.2)+im3;
        
        %% Changing the Intensity Value of Mammogram Image Depend on Histogram of ROI
        
        
        [imm2a, tabela, hha, gga, binLocationsa] = djusting_by_histogram_roi(im3,pathg,pathh,pathbinLocations);
        [imm2b, tabelb, hhb, ggb, binLocationsb] = djusting_by_histogram_roi1(im3,pathg,pathh,pathbinLocations);
        icomba = immultiply(imm2a, (1 - imm2b));
        icomba=imm2a+(imm2b);
        icomba = medfilt2(icomba, [8 8]);
        % Combine images and display
        icomb = immultiply(im4, uint8(icomba));
        imm2 = im4;
        
        
        % Prepare combined image for display
        newImgte = cat(3, im2uint8(im3), im2uint8(im4), im2uint8(icomba));
                figure(1);
                subplot(1, 4, 1, 'align'); imshow(im3, []); xlabel('im3');
                subplot(1, 4, 2, 'align'); imshow(im4, []); xlabel('com1');
                subplot(1, 4, 3, 'align'); imshow(imm2a, []); xlabel('imm2a');
                subplot(1, 4, 4, 'align'); imshow(newImgte, []); xlabel('newImgte');
        
        % Mask processing and image resizing
        [masknewdir] = greatmask11(i1);
        [sz1, sz2] = size(masknewdir);
        [sz11, sz22, sz33] = size(newImgte);
        masknewdir2 = masknewdir;
        
        % Handle left and right images differently
        if imagedir(1, 25) == 'L'  % Left image
            masknewdir2((sz11 + 1):end, :) = [];
            masknewdir2(:, (sz22 + 1):end) = [];
        end
        
        if imagedir(1, 25) == 'R'  % Right image
            masknewdir2((sz11 + 1:sz1), :) = [];
            masknewdir2(:, (sz22 + 1):end) = [];
        end
        
        ZZ = masknewdir2 / 255;
        cond = isequal(size(ZZ), size(newImgte(:,:,1)));
        
        if cond == 1
            Z1 = ZZ .* newImgte;
        else
            disp('The sizes of both images are not the same.');
            continue;
        end
        
        ORIM2 = ~im2double(ZZ) .* im2double(newImgte);
        Z1NEW = removecoulm(Z1);
        
        % Save the processed image
        if typeimage=='B'
            output_filename = sprintf('benignimage%d.jpg', i1);
        end
        
        if typeimage=='N'
            output_filename = sprintf('NORMALimage%d.jpg', i1);
        end
        if typeimage=='C'
            output_filename = sprintf('CANCERimage%d.jpg', i1);
        end
        output_file = fullfile(output_folder, output_filename);
        
        if isempty(Z1NEW)
            disp(['Skipping empty image']);
            continue;  % Skip to the next iteration
        else
            imwrite(Z1NEW, output_file);
        end
        
        
    else
        lop(lo, 1) = i1;
        lo = lo + 1;
    end
end
