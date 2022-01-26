function check = ACR_SliceInversionCheck(img_ACR,obj_ACR)

res_ACR = ACR_RetrievePixelSpacing(obj_ACR);

landmark_1 = (25/rms(res_ACR))/2;
landmark_2 = (85/rms(res_ACR))/2;
[~,r1] = imfindcircles(img_ACR(:,:,1),round([0.8*landmark_1 1.2*landmark_2]),'Sensitivity',0.85);
[~,r11] = imfindcircles(img_ACR(:,:,11),round([0.8*landmark_1 1.2*landmark_2]),'Sensitivity',0.85);

check = 0;
if isempty(r1) && isempty(r11)
    check = 0;
end

if isempty(r1)
    if r11<20
        check = 1;
    else
        check = 0;
    end
elseif isempty(r11)
    if r1 > 20
        check = 1;
    end
elseif r1>20 && r11<20
    check = 1;
else
    check = 0;
end

if check == 0
    centroid = ACR_Centroid(img_ACR,obj_ACR);
    [c1,~] = imfindcircles(img_ACR(:,:,1),round([0.8*landmark_1 1.2*landmark_2]),'Sensitivity',0.95);

    if ~isempty(c1)
        centroid_dist = norm(c1-centroid);
        if centroid_dist < 20
            check = 1;
        else
            check = 0;
        end
    end
end
