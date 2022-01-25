function check = ACR_SliceInversionCheck(img_ACR,obj_ACR)

res_ACR = ACR_RetrievePixelSpacing(obj_ACR);

[~,r1] = imfindcircles(img_ACR(:,:,1),round([10 25]/res_ACR(1)));
[~,r11] = imfindcircles(img_ACR(:,:,11),round([10 25]/res_ACR(1)));

if isempty(r1) && isempty(r11)
    check = 0;
elseif isempty(r11)
    check = 0;
else
    check = 1;
end