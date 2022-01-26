function rot_angle = ACR_FindRotation(img_ACR,obj_ACR)

img = img_ACR(:,:,1);
corr_img = ACR_IntensityGradientCorrection(img,obj_ACR);
b_img = bwareaopen(edge(corr_img,'canny'),600*size(img,1)/600);

[H,theta,~] = hough(b_img,'Theta',-90:0.5:89.5);

P = houghpeaks(H,1);
angle_peak = theta(P(1,2));

if abs(angle_peak) > 45 && abs(abs(angle_peak) - 90) > 0
    if angle_peak < 0 
        rot_angle = 90 + angle_peak;
    else
        rot_angle = -90 + angle_peak;
    end
else
    rot_angle = 0;
end