clc;
clear all;
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');  
addpath('C:\OneDrive - Technion\Phd\MatlabSandbox\Riemmanian Eigenvalues Analysis');
addpath('C:\OneDrive - Technion\Phd\Code\AteamToolbox');
%% Parameters
n = 500;  % 
datafiles_cam1 = 'data/data/s1_%d.jpg'; % Camera 1 image files name format, where %d is the serial number.
datafiles_cam2 = 'data/data/s2_%d.jpg'; % Camera 2 image files name format, where %d is the serial number.
start_at = 100010;     % First image to load
sample_rate = 1;       % Take only one in sample_rate images. 

%% Load images.
sample_rate=1;
disp('Loading')
ResizeFac=0.5;
tic
%
% A sample image, to compute get the size of the vectors. 
%
str0 = sprintf(datafiles_cam1,start_at+1);     % The filename for the snapshot from Camera 1.
Itmp = im2double(imread(str0));     % Read file
% figure(); imshow(Itmp); 
% Itmp=Itmp(:,161:end,:);
Itmp = imresize(Itmp,ResizeFac);        % Downsample

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
    
%     I1=I1(:,161:end,:);
%     I2=I2(:,1:160,:);
%    figure();subplot(1,2,1);imshow(I1); subplot(1,2,2);imshow(I2); 
    
    I1 = imresize(I1,ResizeFac);                % Downsample snapshot 1
    I2 = imresize(I2,ResizeFac);                % Downsample snapshot 2

%    subplot(1,2,1);imshow(I1); subplot(1,2,2);imshow(I2); pause();
    
    ss1(i,:) = reshape(I1, [1, m]);         % Reshape snapshot into a vector, and store it.
    ss2(i,:) = reshape(I2, [1, m]);         % Reshape snapshot into a vector, and store it.
end

%
%   Copy the data array (and modify if needed).
%
s1 = ss1;
s2 = ss2;

toc

%% Distances
disp('Distances')
tic
d1 = pdist2(s1,s1,'euclidean');             % The distance between samples, based only on snapshots taken by Camera 1.     
d2 = pdist2(s2,s2,'euclidean');             % The distance between samples, based only on snapshots taken by Camera 2.  
toc

D1=sqrt(d1);
D2=sqrt(d2);
% figure(); imagesc(d1)
% figure(); imagesc(d2)

%% Parameters
%----data params----
DataParams=struct;
DataParams.NumberOfDatapoints=5e2;

%----kernels construnction params----
NormFac=0.2;
% 0<NormFac<1 - global scale factor that eqauls to NormFac*(median of the pairwise distances)
% NormFac>1 - local scale factor thar equals to the median of the distances  of the NormFac nearest neighbours
UseMutualScaleFac='None';
% UseMutualScaleFac='Mean' - use the same scale factor for each kernel, the mutual scale factor is set to be their mean 
% UseMutualScaleFac='Min' - use the same scale factor for each kernel, the mutual scale factor is set to be the lower scale factor
% UseMutualScaleFac='None' - use different scale factor for each kernel 
KernelsParams=v2struct(NormFac,UseMutualScaleFac);

%----eigenvalues flow diagram params---- 
NumberOfEigenVals=20; %Number of eigenvalues calculated at each point on the geodesic
NumberOfPointAlongTheGeodesicPath=50;% Number of points on the geodesic grid
EvfdParams=v2struct(NumberOfEigenVals,NumberOfPointAlongTheGeodesicPath);

%% Generate datapoints and calculate kernels
v2struct(DataParams);
% [S1,S2,X,Y,Z] = GenerateSyntheticData(DataParams);
S1=s1;S2=s2;
[D1,D2,A1,A2,K1,K2,Scale1,Scale2] = GetKernels(S1,S2,KernelsParams);
% ShowDatapointsAndScaleFactor(S1,S2,X,Y,Z,Scale1,Scale2,KernelsParams);

%% Calculate the eigenvales flow diagram
%----geodeisc interpolator----
[Dim] = EstimateKernelsDim(K1,K2,0);
Interpulator = @ (t) real(FixedGeodesic(K1,K2,t,Dim));

%----linear interpolator----
Interpulator = @ (t) (1-t)*K1+t*K2;

%----harmonic interpolator----
% Interpulator = @ (t) K1^(1-t)*K2^t;

[EigenValuesMat,ColorsMat,tVec]=GetEvfd(EvfdParams,Interpulator);
figure();ShowEvfd(EigenValuesMat,ColorsMat,tVec);

%% Extract angles
str0 = sprintf(datafiles_cam1,start_at+1);     % The filename for the snapshot from Camera 1.
Itmp = im2double(imread(str0));     % Read file
% figure(); imshow(Itmp); 
Itmp=Itmp(:,161:end,:);
Itmp = imresize(Itmp,ResizeFac);        % Downsample

m = size(Itmp, 1)*size(Itmp,2)*size(Itmp,3);

ss1 = zeros(n, m);
ss2 = zeros(n, m);
ss1_str = repmat(str0,n,1);
ss2_str = repmat(str0,n,1);

for i=1:n
    i

    idx = start_at + sample_rate*i;         % Index of sample to be loaded.
    
    str1 = sprintf(datafiles_cam1,idx);     % The filename for the snapshot from Camera 1.
    str2 = sprintf(datafiles_cam2,idx);     % The filename for the snapshot from Camera 2.

    ss1_str(i,:)=str1;                      % Store filename
    ss2_str(i,:)=str2;                      % Store filename
    
    I1 = im2double(imread(str1));           % Load snapshot from Camera 1.
    I2 = im2double(imread(str2));           % Load snapshot from Camera 2.
    
    IYoda=I1(:,1:160,:);
    IBulldog=I1(:,161:end,:);
    IBunny=I2(:,161:end,:);
%    figure();subplot(1,3,1);imshow(IYoda); subplot(1,3,2);imshow(IBulldog);  subplot(1,3,3);imshow(IBunny); 
    
    IYoda = imresize(IYoda,ResizeFac);                % Downsample snapshot 1
    IBulldog = imresize(IBulldog,ResizeFac);                % Downsample snapshot 2
    IBunny = imresize(IBunny,ResizeFac);                % Downsample snapshot 2
    
    ss_yoda(i,:) = reshape(IYoda, [1, m]);         % Reshape snapshot into a vector, and store it.
    ss_bulldog(i,:) = reshape(IBulldog, [1, m]);         % Reshape snapshot into a vector, and store it.
    ss_bunny(i,:) = reshape(IBunny, [1, m]);         % Reshape snapshot into a vector, and store it.

end

d_yoda = pdist2(ss_yoda,ss_yoda,'euclidean');             % The distance between samples, based only on snapshots taken by Camera 1.     
d_bulldog = pdist2(ss_bulldog,ss_bulldog,'euclidean');             % The distance between samples, based only on snapshots taken by Camera 2.  
d_bunny = pdist2(ss_bunny,ss_bunny,'euclidean');             % The distance between samples, based only on snapshots taken by Camera 2.  

% figure(); imagesc(d1)
% figure(); imagesc(d2)
Dyoda=sqrt(d_yoda);
Dbulldog=sqrt(d_bulldog);
Dbunny=sqrt(d_bunny);

% NormFac=3;
epsyoda=median(Dyoda(:))*NormFac;
epsbulldog=median(Dbulldog(:))*NormFac;
epsbunny=median(Dbunny(:))*NormFac;


epsyoda=4.5;
epsbulldog=4.5;
epsbunny=4.5;

Ayoda=exp(-Dyoda.^2/(epsyoda^2));
Kyoda = sinkhornKnopp(Ayoda);Kyoda=Kyoda^1;
Abulldog=exp(-Dbulldog.^2/(epsbulldog^2));
Kbulldog = sinkhornKnopp(Abulldog);Kbulldog=Kbulldog^1;
Abunny=exp(-Dbunny.^2/(epsbunny^2));
Kbunny = sinkhornKnopp(Abunny);Kbunny=Kbunny^1;

[MapEmbdbunny, ~, ~,singvalsbunny, ~,~] =DiffusionMapsFromKer( Kbunny , 1);
[MapEmbdbulldog, ~, ~,singvalsbulldog, ~,~] =DiffusionMapsFromKer( Kbulldog , 1);
[MapEmbdyoda, ~, ~,singvalsyoda, ~,~] =DiffusionMapsFromKer( Kyoda , 1);

numNeighbors_refine=20;
tau=1;
alg1_d_t = pdist2(Kbunny',Kbunny','euclidean');
[MapEmbdbunny, UWeighted, d_tau, alg1_singvals]=DiffusionMapsFromDistance(alg1_d_t,tau,numNeighbors_refine);
alg1_d_t = pdist2(Kbulldog',Kbulldog','euclidean');
[MapEmbdbulldog, UWeighted, d_tau, alg1_singvals]=DiffusionMapsFromDistance(alg1_d_t,tau,numNeighbors_refine);
alg1_d_t = pdist2(Kyoda',Kyoda','euclidean');
[MapEmbdyoda, UWeighted, d_tau, alg1_singvals]=DiffusionMapsFromDistance(alg1_d_t,tau,numNeighbors_refine);


% figure(); scatter(MapEmbdyoda(:,2),MapEmbdyoda(:,3))
% figure(); scatter(MapEmbdbulldog(:,2),MapEmbdbulldog(:,3))
% figure(); scatter(MapEmbdbunny(:,2),MapEmbdbunny(:,3))
% 
% figure(); scatter3(MapEmbdyoda(:,2),MapEmbdyoda(:,3),MapEmbdyoda(:,4))
% figure(); scatter3(MapEmbdbulldog(:,2),MapEmbdbulldog(:,3),MapEmbdbulldog(:,4))
% figure(); scatter3(MapEmbdbunny(:,2),MapEmbdbunny(:,3),MapEmbdbunny(:,4))
% 

AngleYoda=atan(abs(MapEmbdyoda(:,2)./MapEmbdyoda(:,3)));
% figure(); scatter(MapEmbdyoda(:,2),MapEmbdyoda(:,3),10,AngleYoda)
AngleBulldog=atan(abs(MapEmbdbulldog(:,2)./MapEmbdbulldog(:,3)));
AngleBunny=atan(abs(MapEmbdbunny(:,2)./MapEmbdbunny(:,3)));


%%
GeodesicEigenValuesRawStack=log(abs(EigenValuesMat(:)));
tRawStack=kron(tVec(1:size(EigenValuesMat,2)),ones(1,size(EigenValuesMat,1)))';
GeodesicEigenValuesTupple=[GeodesicEigenValuesRawStack,tRawStack];


VideoNameFile=["SweepSimulationPaper"];
v = VideoWriter(char(VideoNameFile+".avi"));
v.FrameRate=5;
open(v);

figure();
set(gcf, 'Position', get(0, 'Screensize'));

for t0=tVec
    subplot(3,6,[1:3,7:9,13:15])
    ShowEvfd(EigenValuesMat,ColorsMat,tVec);
    title('Eigenvalues flow diagram','FontSize',30);
    ylabel('$t$','FontSize',30);xlabel('$-\log(\mu_t^i)$','FontSize',30);
    [~,RelevantRow]=min(abs(t0-tVec));
    hold on;
    leadingevs=(NumberOfEigenVals-1)*(RelevantRow-1)+(1:3);
    scatter((GeodesicEigenValuesTupple(leadingevs,1)),GeodesicEigenValuesTupple(leadingevs,2),200,'r','o','LineWidth',4);
    hold on;
    plot([min(GeodesicEigenValuesTupple(:,1)),max(GeodesicEigenValuesTupple(:,1))],[tVec(RelevantRow),tVec(RelevantRow)],'r','LineWidth',2)
    axis tight;
%     pause()
    hold off;

    K=FixedGeodes( K1,K2,tVec(RelevantRow),-1 );K=real(K);
    [MapEmbd, ~, ~,singvals, ~,~] =DiffusionMapsFromKer( K , 1);
    
    subplot(3,6,3+1);
    scatter(MapEmbd(:,2),MapEmbdyoda(:,2),'filled');
%     xlabel('$\Psi_1$');
    ylabel('$v_t^2$','FontSize',30);set(gca,'color','none')
    subplot(3,6,3+2);
    scatter(MapEmbd(:,2),MapEmbdbulldog(:,2),'filled');
%     xlabel('$\Psi_1$');ylabel('$\Psi_1^{Bulldog}$');
    set(gca,'color','none')
    subplot(3,6,3+3);
    scatter(MapEmbd(:,2),MapEmbdbunny(:,2),'filled');
%     xlabel('$\Psi_1$');ylabel('$\Psi_1^{Bunny}$');
    set(gca,'color','none')
    %
    subplot(3,6,6+4);
    scatter(MapEmbd(:,3),MapEmbdyoda(:,2),'filled');
%     xlabel('$\Psi_2$');ylabel('$\Psi_1^{Yoda}$');
    ylabel('$v_t^3$','FontSize',30);set(gca,'color','none')
    subplot(3,6,6+5);
    scatter(MapEmbd(:,3),MapEmbdbulldog(:,2),'filled');
%     xlabel('$\Psi_2$');ylabel('$\Psi_1^{Bulldog}$');
    set(gca,'color','none')
    subplot(3,6,6+6);
    scatter(MapEmbd(:,3),MapEmbdbunny(:,2),'filled');
%     xlabel('$\Psi_2$');ylabel('$\Psi_1^{Bunny}$');
    set(gca,'color','none')
    %
    subplot(3,6,9+7);
    scatter(MapEmbd(:,4),MapEmbdyoda(:,2),'filled');
%     xlabel('$\Psi_3$');ylabel('$\Psi_1^{Yoda}$');
    ylabel('$v_t^4$','FontSize',30);set(gca,'color','none')
    xlabel('$\Theta_{yoda}$','FontSize',30);
    subplot(3,6,9+8);
    scatter(MapEmbd(:,4),MapEmbdbulldog(:,2),'filled');
%     xlabel('$\Psi_3$');ylabel('$\Psi_1^{Bulldog}$');
    xlabel('$\Theta_{bulldog}$','FontSize',30);
    set(gca,'color','none')
    subplot(3,6,9+9);
    scatter(MapEmbd(:,4),MapEmbdbunny(:,2),'filled');
%     xlabel('$\Psi_3$');ylabel('$\Psi_1^{Bunny}$');
    xlabel('$\Theta_{bunny}$','FontSize',30);
    set(gca,'color','none')
    
    axes('pos',[.55 .92 .08 .08]);imshow(IYoda)
    axes('pos',[.67 .92 .08 .08]);imshow(IBulldog)
    axes('pos',[.81 .92 .08 .08]);imshow(IBunny)
    writeVideo(v,getframe(gcf));
end
close(gcf);
close(v);

