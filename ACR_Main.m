%% ACR Automation Main Script
% BY YASSINE AZMA (OCT 2021)
%
% yassine.azma@rmh.nhs.uk
% yassine.azma01@icr.ac.uk
%
% If preferable, please raise issues on 
% github.com/YassineRMH/MATLAB_ACR_Automation
%
% All scripts are currently configured to run with both single-frame
% and multi-frame data. It has been tested for Siemens, GE and Philips
% datasets. It is not known if these scripts work for other manufacturers.
%
% This suite of scripts require the following toolboxes:
% 
% Image Processing Toolbox
% Signal Processing Toolbox
% Curve Fitting Toolbox ** only required for MTF
%
%
%%%%%%%%%%                      OPTIONS                          %%%%%%%%%%
%%% Physical Orientation of Phantom - CURRENTLY ONLY NEEDED FOR LOCALISER
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

%%% Rician Noise Correction
% Selection ['yes' or 'no']
%
% Opt to include a correction factor to single-image SNR values for Rician
% noise.
options.RicianNoiseCorrection = 'yes';

%%% Slice Displacement for Slice Position
% Selection ['yes' or 'no']
%
% Opt to report slice displacement in the workspace instead of the standard 
% bar offset. For reference, this is considered by the ACR guidance as half
% the bar offset for each respective slice.
options.SliceDisplacementDisplay = 'no';

%%% Directory for Spreadsheet Output
% Selection ['PATH']
% 
% Select the directory where the spreadsheet containing results will be 
% saved. If left blank, the spreadsheet will default to the current active 
% directory in the Current Folder section of the UI.
options.OutputPath = 'D:\STP\MRI\ACR Data\ACR Results'

%% Data Sorting
% Upon running this script, you will be prompted to open data from a folder 
% containing your sagittal localiser. If your localiser folder contains 
% multiple images (i.e a multiplanar localiser), the first sagittal image 
% will be chosen. 
% 
% You will also be asked to provide either one or two axial series.
% Select 'Two' if you would like to perform the NEMA subtraction method for
% calculating SNR.

clearvars -except options
[img_loc,img_ACR,obj_loc,obj_ACR] = ACR_DataSort(options);
%% Geometric Accuracy

L = ACR_GeometricAccuracy(img_loc,img_ACR,obj_loc,obj_ACR)
%% High-Contrast Spatial Resolution

resolvable = ACR_SpatialResolution(img_ACR,obj_ACR)
%% Slice Thickness Accuracy

dz = ACR_SliceThickness(img_ACR,obj_ACR)
%% Slice Position Accuracy

dL = ACR_SlicePosition(img_ACR,obj_ACR,options) % slice 1 and 11
%% Image Intensity Uniformity
    
PIU = ACR_Uniformity(img_ACR,obj_ACR)
%% Percent-Signal Ghosting

PSG = ACR_Ghosting(img_ACR,obj_ACR)
%% Low-Contrast Object Detectability

% DO IT MANUALLY
%% SNR

SNR = ACR_SNR(img_ACR,obj_ACR,options)
%% SNR NEMA Subtraction

% if size(img_ACR,4) > 1
%     sub_SNR = ACR_subSNR(img_ACR,obj_ACR)
% else
%     sub_SNR = [];
% end
%% MTF (Experimental) ONLY USE ON SLANTED EDGE

% eff_res = ACR_MTF(img_ACR,obj_ACR)

%% Report

ACR_Report(L,resolvable,dz,dL,PIU,PSG,SNR,obj_ACR,obj_loc,options)