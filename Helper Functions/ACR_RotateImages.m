function rot_img = ACR_RotateImages(img_ACR,rot_angle)

if rot_angle ~= 0
    disp(['Phantom tilt of ' num2str(rot_angle) char(176) ' detected. Rotating images...'])
    for k = 1:size(img_ACR,3)
        ss_img = imresize(img_ACR(:,:,k),4); % supersample
        rot_ss_img = imrotate(ss_img,rot_angle,'crop','bicubic'); % perform rotation
        rot_img(:,:,k) = imresize(rot_ss_img,1/4); % downsample
    end
else
    rot_img = img_ACR;
end