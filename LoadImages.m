%% Unpack and load images
datafiles_cam1 = 'Images/s1_%d.jpg'; % Camera 1 image files name format, where %d is the serial number.
datafiles_cam2 = 'Images/s2_%d.jpg'; % Camera 2 image files name format, where %d is the serial number.

n = DataParams.NumberOfDatapoints;  %
start_at = DataParams.start_at;     %
sample_rate = DataParams.sample_rate;    
ResizeFac=DataParams.ResizeFac;



tic
%
% A sample image, to compute get the size of the vectors. 
%
str0 = sprintf(datafiles_cam1,start_at+1);     % The filename for the snapshot from Camera 1.
Itmp = im2double(imread(str0));                % Read file
Itmp = imresize(Itmp,ResizeFac);               % Downsample

%
% Image size
%
m = size(Itmp, 1)*size(Itmp,2)*size(Itmp,3);

%
% Initialize arrays.
% ss1 and ss2 are going to store the downsampled images, 
% each image reshaped into a vector.
% ss1_str and ss2_str are going to store the filenames from which we read
% the images.
%
ss1 = zeros(n, m);
ss2 = zeros(n, m);
ss1_str = repmat(str0,n,1);
ss2_str = repmat(str0,n,1);

%
% Load images one by one
%
for i=1:n
    %i

    idx = start_at + sample_rate*i;         % Index of sample to be loaded.
    
    str1 = sprintf(datafiles_cam1,idx);     % The filename for the snapshot from Camera 1.
    str2 = sprintf(datafiles_cam2,idx);     % The filename for the snapshot from Camera 2.

    ss1_str(i,:)=str1;                      % Store filename
    ss2_str(i,:)=str2;                      % Store filename
    
    I1 = im2double(imread(str1));           % Load snapshot from Camera 1.
    I2 = im2double(imread(str2));           % Load snapshot from Camera 2.
        
    I1 = imresize(I1,ResizeFac);                % Downsample snapshot 1
    I2 = imresize(I2,ResizeFac);                % Downsample snapshot 2
    
    ss1(i,:) = reshape(I1, [1, m]);         % Reshape snapshot into a vector, and store it.
    ss2(i,:) = reshape(I2, [1, m]);         % Reshape snapshot into a vector, and store it.
end

%
%   Copy the data array (and modify if needed).
%
s1 = ss1;
s2 = ss2;

toc
