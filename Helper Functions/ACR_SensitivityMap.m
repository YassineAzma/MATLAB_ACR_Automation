%% ACR Sensitivity Map
% by Yassine Azma (Apr 2022)
%
% This script creates a rough sensitivity map based on slice 7 of the ACR
% phantom.

function sens_map = ACR_SensitivityMap(img_ACR)

temp = img_ACR(:,:,7);
max_val = max(temp(:));

sens_map = temp./max_val;
sens_map = imguidedfilter(sens_map,'DegreeOfSmoothing',1);