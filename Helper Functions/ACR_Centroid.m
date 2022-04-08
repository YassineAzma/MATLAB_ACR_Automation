%% ACR Data Sort
% by Yassine Azma (Jan 2022)
% 
% This script calculates the centroid of the ACR phantom using a 
% centre-of-mass calculation.

function centroid = ACR_Centroid(img_ACR,obj_ACR)

igc_img_ACR = ACR_IntensityGradientCorrection(img_ACR,obj_ACR);

img_insert = igc_img_ACR(:,:,7);
thresh = bwareaopen(img_insert>0.25*max(img_insert(:)),500); % threshold and remove unconnected pixels
bhull = bwconvhull(thresh); % create convex hull image

% Centre of Mass Method
[row,col] = find(bhull);
centroid = round([sum(col)/length(col), sum(row)/length(row)]);

% Circular Arc Method (for incomplete thresholding)
% arc = bwareaopen(edge(bhull,'canny',0.5),500);
% [row,col] = find(arc);
% coord_try = [randi(length(row),5,1)];
% 
% coord_perms = nchoosek(coord_try,3);
% 
% for k = 1:size(coord_perms,1)
%     ax = col(coord_perms(k,1));
%     bx = col(coord_perms(k,2));
%     cx = col(coord_perms(k,3));
% 
%     ay = row(coord_perms(k,1));
%     by = row(coord_perms(k,2)); 
%     cy = row(coord_perms(k,3));
% 
%     D = [ax-bx, ay-by;
%         bx-cx, by-cy];
%     T = 0.5*[ax^2-bx^2 + ay^2-by^2
%         bx^2-cx^2 + by^2-cy^2];
% 
%     c(k,:) = D\T;
% end
% 
% centroid_arc = round(mean(c));