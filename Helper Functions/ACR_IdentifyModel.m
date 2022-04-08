%% ACR Model Identification
% by Yassine Azma (Apr 2022)
%
% This script identifies whether the ACR phantom is the newer model with a 
% rectangular void board vs the older model with the distortion grid.

function model = ACR_IdentifyModel(img_ACR)

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_insert = squeeze(double(img_ACR(:,:,5,1))); % if yes, only process the first
    waitfor(msgbox('4D array detected. Only processing first axial series.'));
else
    img_insert = double(img_ACR(:,:,5));
end

thresh_img = img_insert > 0.2*max(img_insert(:));
label_map = bwlabel(thresh_img,4);

if max(label_map(:)) > 15
    model = 'old';
else
    model = 'new';
end