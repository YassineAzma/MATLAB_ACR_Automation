%% ACR Ghosting
% by Yassine Azma (Oct 2021)
%

function PSG = ACR_Ghosting(img_ACR,obj_ACR,options)

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_unif = double(img_ACR(:,:,7,:)); % if yes, only process the first
else
    img_unif = double(img_ACR(:,:,7));
end

for k = 1:size(img_ACR,4)
    res_ACR = ACR_RetrievePixelSpacing(obj_ACR(k));

    r = 8; % for ~201cm^2 circular ROI
    r_img = ceil(10*r./res_ACR(1)); % equivalent pixel radius
    d_void = ceil(5/res_ACR(1));

    % Find centroid
    centroid = ACR_Centroid(img_ACR(:,:,:,k),obj_ACR(k)); % determine centroid from convex hull image

    mask = ACR_Threshold(img_ACR(:,:,:,k),res_ACR,centroid);
    temp = img_unif(:,:,k);

    % Large ROI
    [img_cols, img_rows] = meshgrid(1:size(temp,1),1:size(temp,2)); % create grid
    roi_index = (img_rows - centroid(2) - d_void).^2 + (img_cols - centroid(1)).^2 <= r_img.^2; % create large ROI

    % Elliptical ROIs
    w_point = find(sum(mask,1)>0,1,'first')-1; % westmost point
    if w_point > 0.1*size(temp,1) || w_point > centroid(1)
        w_point = centroid(1)-95/res_ACR(1);
    end
    e_point = find(sum(mask,1)>0,1,'last')-1; % eastmost point
    if e_point < 0.9*size(temp,1) || e_point < centroid(1)
        e_point = centroid(1)+95/res_ACR(1);
    end
    n_point = find(sum(mask,2)>0,1,'first')-1; % northmost point
    if n_point > 0.1*size(temp,2) || n_point > centroid(2)
        n_point = centroid(2)-95/res_ACR(1);
    end
    s_point = find(sum(mask,2)>0,1,'last')-1; % southmost point
    if s_point < 0.9*size(temp,2) || s_point < centroid(2)
        s_point = centroid(2)+95/res_ACR(1);
    end

    w_ellip_centre = [centroid(2) floor(w_point/2)];
    w_ellip_roi_index = ((img_rows - w_ellip_centre(1))/4).^2 + (img_cols - w_ellip_centre(2)).^2 <= ceil(10./res_ACR(1)).^2; % create W large ROI

    e_ellip_centre = [centroid(2) (e_point + ceil((size(img_unif,2)-e_point)/2))];
    e_ellip_roi_index = ((img_rows - e_ellip_centre(1))/4).^2 + (img_cols - e_ellip_centre(2)).^2 <= ceil(10./res_ACR(1)).^2; % create E large ROI

    n_ellip_centre = [round(n_point/2) centroid(1)];
    n_ellip_roi_index = ((img_rows - n_ellip_centre(1))).^2 + ((img_cols - n_ellip_centre(2))/4).^2 <= ceil(10./res_ACR(1)).^2; % create N large ROI

    s_ellip_centre = [(s_point + round((size(img_unif,2)-s_point)/2)) centroid(1)];
    s_ellip_roi_index = (img_rows - s_ellip_centre(1)).^2 + ((img_cols - s_ellip_centre(2))/4).^2 <= ceil(10./res_ACR(1)).^2; % create S large ROI

    large_roi_val = mean(nonzeros(temp(roi_index)));
    w_ellip_roi_val = mean(nonzeros(temp(w_ellip_roi_index)));
    e_ellip_roi_val = mean(nonzeros(temp(e_ellip_roi_index)));
    n_ellip_roi_val = mean(nonzeros(temp(n_ellip_roi_index)));
    s_ellip_roi_val = mean(nonzeros(temp(s_ellip_roi_index)));

    PSG(k) = 100*abs(((n_ellip_roi_val + s_ellip_roi_val)-(w_ellip_roi_val + e_ellip_roi_val))/(2*large_roi_val));

    if strcmp(options.SuppressFigures,'no')
        figure
        imshow(temp,[])
        hold on
        plot(r_img*cosd(0:1:360)+centroid(1),r_img*sind(0:1:360)+centroid(2)+5)
        text(centroid(1)-3*floor(10./res_ACR(1)), centroid(2)+floor(10./res_ACR(1)),...
            sprintf('mean = %.1f',large_roi_val),'color','w','fontsize',10) % label with measured mean

        plot(ceil(10./res_ACR(1))*cosd(0:1:360)+w_ellip_centre(2),ceil(10./res_ACR(1))*4*sind(0:1:360)+w_ellip_centre(1),'color','r')
        text(w_ellip_centre(2), w_ellip_centre(1)+floor(10./res_ACR(1)),...
            sprintf('mean = %.1f',w_ellip_roi_val),'color','w','fontsize',10) % label with measured mean

        plot(ceil(10./res_ACR(1))*4*cosd(0:1:360)+n_ellip_centre(2),ceil(10./res_ACR(1))*sind(0:1:360)+n_ellip_centre(1),'color','r')
        text(n_ellip_centre(2)-5*floor(10./res_ACR(1)), n_ellip_centre(1),...
            sprintf('mean = %.1f',n_ellip_roi_val),'color','w','fontsize',10) % label with measured mean

        plot(ceil(10./res_ACR(1))*cosd(0:1:360)+e_ellip_centre(2),ceil(10./res_ACR(1))*4*sind(0:1:360)+e_ellip_centre(1),'color','r')
        text(e_ellip_centre(2)-5*floor(10./res_ACR(1)), e_ellip_centre(1),...
            sprintf('mean = %.1f',e_ellip_roi_val),'color','w','fontsize',10) % label with measured mean

        plot(ceil(10./res_ACR(1))*4*cosd(0:1:360)+s_ellip_centre(2),ceil(10./res_ACR(1))*sind(0:1:360)+s_ellip_centre(1),'color','r')
        text(s_ellip_centre(2), s_ellip_centre(1),...
            sprintf('mean = %.1f',s_ellip_roi_val),'color','w','fontsize',10) % label with measured mean
        % xlabel(['PSG = ' num2str(round(PSG,3)) '%'],'fontweight','bold','fontsize',14)
        % title('Ghosting')
    end
end

if isrow(PSG)
    PSG = PSG';
end