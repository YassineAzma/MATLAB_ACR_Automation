%% ACR Automation Main Script
% It is highly recommended to run the script one section at a time using
% Ctrl+Enter. All scripts are currently configured to run with both single-frame
% and multi-frame data.
%
% This suite of scripts require the following toolboxes:
% 
% Image Processing Toolbox
% Signal Processing Toolbox

%% Data Sorting
% Upon running this script, you will be prompted to open data from a folder 
% containing your sagittal localiser. If your localiser folder contains 
% multiple images (i.e a multiplanar localiser), the first sagittal image 
% will be chosen. 
% 
% You will also be asked to provide either one or two axial series.
% Select 'Two' if you would like to perform the NEMA subtraction method for
% calculating SNR.

clearvars

% OPTIONS
%%% Physical Orientation of Phantom 
% Selection: ['axial','sagittal' or 'coronal']
%
% This should be 'axial' unless the phantom is set up in an orientation 
% outside of the ACR guidance.
options.Orientation = 'axial'; 

%%% Include Localiser 
% Selection: ['yes' or 'no']
%
% Opt to include the localiser in the geometric accuracy test
options.IncludeLocaliser = 'no';

[img_loc,img_ACR,obj_loc,obj_ACR,info] = ACR_DataSort(options);
%% Geometric Accuracy

L = ACR_GeometricAccuracy(img_loc,img_ACR,obj_loc,obj_ACR)
%% High-Contrast Spatial Resolution

resolvable = ACR_SpatialResolution(img_ACR,obj_ACR)
%% Slice Thickness Accuracy

dz = ACR_SliceThickness(img_ACR,obj_ACR)
%% Slice Position Accuracy

dL = ACR_SlicePosition(img_ACR,obj_ACR) % slice 1 and 11
%% Image Intensity Uniformity
    
PIU = ACR_Uniformity(img_ACR,obj_ACR)
%% Percent-Signal Ghosting

PSG = ACR_Ghosting(img_ACR,obj_ACR)
%% Low-Contrast Object Detectability
% DO IT MANUALLY
%% SNR

SNR = ACR_SNR(img_ACR,obj_ACR)
%% SNR NEMA Subtraction

if size(img_ACR,4) > 1
    sub_SNR = ACR_subSNR(img_ACR,obj_ACR)
end