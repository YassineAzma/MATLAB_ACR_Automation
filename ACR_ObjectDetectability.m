%% ACR Low Contrast Object Detectability
% by Yassine Azma (Nov 2021)
% 
% This script ...

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_spokes = squeeze(double(img_ACR(:,:,9:11,1))); % if yes, only process the first
    waitfor(msgbox('4D array detected. Only processing first axial series.'));
else
    img_spokes = double(img_ACR(:,:,9:11));
end

res_ACR = ACR_RetrievePixelSpacing(obj_ACR);

centroid_mask = img_spokes(:,:,1) > 0.15*max(img_spokes(:)); % convex hull image
centroid = regionprops(centroid_mask,'Centroid'); % find centroid
centroid = centroid(2).Centroid; % take inner circle centroid
spoke_rad = round(45.3/rms(res_ACR)); % in pixels

[img_cols, img_rows] = meshgrid(1:size(img_spokes,1),1:size(img_spokes,2)); % create grid 
roi_index = (img_rows - centroid(2)).^2 + (img_cols - centroid(1)).^2 <= 0.2*spoke_rad.^2; % create large ROI 

init_angle = [1:180];
temp_img1 = img_spokes(:,:,1);
temp_img2 = img_spokes(:,:,2);
temp_img3 = img_spokes(:,:,3);
baseline = [min(nonzeros(temp_img1(roi_index))) min(nonzeros(temp_img2(roi_index))) min(nonzeros(temp_img3(roi_index)))];

