clc;
clear all;
addpath('Utils\')

%% Load parameters
Parameters;

%% Load data
disp('Unzipping  images folder')
unzip('Images.zip')

disp('Loading images from folder')
LoadImages;

%% Generate datapoints and calculate kernels
v2struct(DataParams);

[D1,D2,A1,A2,K1,K2,Scale1,Scale2] = GetKernels(s1,s2,KernelsParams);

%% Calculate the eigenvales flow diagram
v2struct(EvfdParams);

switch Interpolator
    case 'Geodesic'
        [Dim] = EstimateKernelsDim(K1,K2,TolFac);
        Interpulator = @ (t) FixedGeodesic(K1,K2,t,Dim);
    case 'Linear'
        Interpulator = @ (t) (1-t)*K1+t*K2;
    case 'Harmonic'
        Interpulator = @ (t) K1^(1-t)*K2^t;
end

[EigenValuesMat,ColorsMat,tVec]=GetEvfd(EvfdParams,Interpulator);
figure();ShowEvfd(EigenValuesMat,ColorsMat,tVec);


