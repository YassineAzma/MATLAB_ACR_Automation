%% ACR Intensity Gradient Correction
% by Yassine Azma (Dec 2021)
%
% This script takes the ACR data and models it as a flat constant value
% plus added orthogonal gradients and parabolas. Solving the regression 
% problem using weighted linear least squares allows for the production of 
% an image which is much flatter. This is used in the ACR_Centroid script.

function corr_img_ACR = ACR_IntensityGradientCorrection(img_ACR)

dims = size(img_ACR,[1 2]);

I0 = ones(dims); % offset
Ix = repmat(-0.5:1/(dims(1)-1):0.5,dims(2),1); % col gradient
Ix2 = Ix.*Ix; % quadratic 
Iy = repmat(-0.5:1/(dims(1)-1):0.5,dims(2),1)'; % row gradient
Iy2 = Iy.*Iy; % quadratic

H = [I0(:) Ix(:) Iy(:) Ix2(:) Iy2(:)]; % model matrix

for k = 1:size(img_ACR,3)
    img = img_ACR(:,:,k);
    wm = im2double(img)/im2double(max(img,[],'all')); % weighting based on mag
    wv = wm(:); % vectorise

    HpWH = H'*((wv*ones(1,5)).*H); % First order
    aw = HpWH \ (H'*(wv.*img(:))); % Ordinary least squares for coefficients

    if abs(aw(2)) > abs(aw(3))
        corr_img = img-imfill(wm>0.03,'holes').*reshape(H(:,[2 4])*aw([2,4]),dims);
        edge_img = bwconvhull(bwareaopen(edge(corr_img,'Canny'),500));
        corr_img_ACR(:,:,k) = edge_img.*corr_img;
    else
        corr_img = img-imfill(wm>0.03,'holes').*reshape(H(:,[3 5])*aw([3,5]),dims);
        edge_img = bwconvhull(bwareaopen(edge(corr_img,'Canny'),500));
        corr_img_ACR(:,:,k) = edge_img.*corr_img;
    end
end