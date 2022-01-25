%% ACR Threshold
% by Yassine Azma (Jan 2022)
%
% This script takes the (axial) ACR series and uses an edge detection
% filter in conjunction with binary labelling to identify the pixel group
% that coincides with the contour of the phantom. The correct pixel group
% is identified using the known radius of the phantom (and the centroid)
% and a binary convex hull image is then created as an image
% mask.

function mask = ACR_Threshold(img,res,centroid)

for n = 1:size(img,3)
    temp = img(:,:,n);
    edge_img = bwareaopen(edge(temp,'Canny'),1500*size(img,1)/600); 
    label_map = bwlabel(edge_img);
    if max(label_map,[],'all') > 1
        for k = 1:max(label_map,[],'all')
            label_img = label_map == k;
            label_centroid(:,k) = regionprops(label_img).Centroid;
            label_bbox(:,k) = regionprops(label_img).BoundingBox;
            label_radius(k) = rms(res)*mean(label_bbox(3:4,k))/2;

            if nargin == 3
                centroid_test(k) = rms(label_centroid(:,k) - centroid');
            end
            radius_test(k) = abs(label_radius(k) - 95);
        end
        if nargin == 3
            test = sqrt(centroid_test.^2+radius_test.^2);
        else
            test = radius_test;
        end
        ind = find(test==min(test));
        mask = bwconvhull(label_map==ind);
    else
        mask = bwconvhull(edge_img);
    end
end