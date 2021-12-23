%% ACR Orientation Check
% by Yassine Azma (Dec 2021)
%
% This script takes sagittal or coronal data (in place of the axial series) 
% and compares it to an exemplar axial set in order to perform necessary
% rotations and slice order corrections. This is only useful for centres
% who perform the ACR tests in multiple orientations, rather than only
% axially as recommended in the guidance.

function [img_ACR,info] = ACR_OrientationCheck(img_ACR,obj_ACR,options)

exemp = load('exemplar_axial.mat').exemp;

if isempty(obj_ACR.getAttributeByName('PixelSpacing')) % Multi-frame check
    list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');
    res_ACR = list.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
else
    res_ACR = obj_ACR.getAttributeByName('PixelSpacing'); % retrieve ACR in-plane resolution
end

%% Rotations

% Intensity Correction

if strcmp(options.IntensityCorrection,'yes') 
    corr_img_ACR = ACR_IntensityGradientCorrection(img_ACR);
else
    corr_img_ACR = img_ACR;
end

mask = bwconvhull(bwareaopen(corr_img_ACR(:,:,6) > 0.2*max(corr_img_ACR(:,:,6),[],'all'),50)); % threshold and make convex hull
centroid = regionprops(mask,'Centroid').Centroid; % determine centroid from convex hull image

w_point = find(sum(mask,1)>0,1,'first')-1; % westmost point
e_point = find(sum(mask,1)>0,1,'last')-1; % eastmost point
n_point = find(sum(mask,2)>0,1,'first')-1; % northmost point
s_point = find(sum(mask,2)>0,1,'last')-1; % southmost point

line_prof_v = improfile(corr_img_ACR(:,:,6),[centroid(1) centroid(1)],[n_point s_point]); % take a vertical line profile
line_prof_h = improfile(corr_img_ACR(:,:,6),[w_point e_point],[centroid(2) centroid(2)]); % take a horizontal line profile

[~,locs_v] = findpeaks(line_prof_v,'MinPeakDistance',25/res_ACR(2),'MinPeakProminence',0.5*max(line_prof_v)); % top to bottom
[~,locs_h] = findpeaks(line_prof_h,'MinPeakDistance',25/res_ACR(1),'MinPeakProminence',0.5*max(line_prof_h)); % left to right
[~,locs_vr] = findpeaks(flip(line_prof_v),'MinPeakDistance',25/res_ACR(2),'MinPeakProminence',0.5*max(line_prof_v)); % bottom to top
[~,locs_hr] = findpeaks(flip(line_prof_h),'MinPeakDistance',25/res_ACR(1),'MinPeakProminence',0.5*max(line_prof_h)); % right to left

min_list = [min(locs_v),min(locs_h),min(locs_vr),min(locs_hr)];
indices = find(min_list); % find minimum distance to first peak as representative of slice position insert
ind = find(min_list == min(min_list(indices)));

switch ind
    case 2
        r_img_ACR = imrotate(img_ACR,-90,'bilinear','crop');
        corr_img_ACR = imrotate(corr_img_ACR,-90,'bilinear','crop');
    case 3
        r_img_ACR = imrotate(img_ACR,-180,'bilinear','crop');
        corr_img_ACR = imrotate(corr_img_ACR,-180,'bilinear','crop');
    case 4
        r_img_ACR = imrotate(img_ACR,-270,'bilinear','crop');
        corr_img_ACR = imrotate(corr_img_ACR,-270,'bilinear','crop');
end

%% Image Registration

% Registration
[optimizer,metric]=imregconfig('multimodal');

reg_img_ACR_1 = imregister(corr_img_ACR(:,:,1),exemp(:,:,1),'affine',optimizer,metric);
reg_img_ACR_2 = imregister(corr_img_ACR(:,:,11),exemp(:,:,1),'affine',optimizer,metric);

%% Similarity
sim = [ssim(reg_img_ACR_1,exemp(:,:,1)),ssim(reg_img_ACR_2,exemp(:,:,1))];

% Reorder

switch ind
    case 2
        img_ACR = imrotate(img_ACR,-90,'bilinear','crop');
    case 3
        img_ACR = imrotate(img_ACR,-180,'bilinear','crop');
    case 4
        img_ACR = imrotate(img_ACR,-270,'bilinear','crop');
end

if sim(2) > sim(1)
    img_ACR = flip(r_img_ACR,3);
end

info.centroid = centroid;