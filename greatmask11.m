function [Z2 ] = greatmask11(i )
 
%%
%% read image

% %%imshow( (maskim))BW2 = imfill(maskim,'holes');

%%
x=0;
h=double(zeros(256,1));
g=double(zeros(256,1));
[path1]=get_excel_file_paths()
A=readtable(path1);
ihi=0;
iihi=0;

imagedir = char(fullfile('D:\',A.fullPath(i)));
cancerimage=imread(imagedir);
maskdir = char(fullfile('D:\',A.Tumour_Contour(i)));
maskim=imread(maskdir);
% BW2 = imfill(maskim,'holes');

%%
[sz1, sz2]=size(cancerimage);
[sz3, sz4]=size(maskim);
zeroim=(zeros(sz1,sz2));
cal1=sz1*sz2;
cal2=sz3*sz4;
if cal1>cal2
    zeroim(1:sz3,1:sz4)=maskim;
    
elseif cal1<cal2
    zeroim=maskim(1:sz1,1:sz2);
    
else
    zeroim=maskim;
end


BW2 = imfill(zeroim,'holes');

ZZ=255-BW2;


ZZ=ZZ/255;
Z1=ZZ.*cancerimage;
Z2=Z1;
Z2=255-(Z2*255);

Z1(1:100,:) = [];
Z1(2100:end,:) = [];

Z1=RemoveLabColum(Z1);
if isempty(Z1)
    
    img=imbinarize(Z1,0.1);
 
else
img=Z1*255;
img=255-img;

end
outputArg1  = img;
end

