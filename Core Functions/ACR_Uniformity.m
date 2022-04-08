%% ACR Uniformity
% by Yassine Azma (Oct 2021)
%
% This script uses the ACR series to place small ROIs (~1cm^2) within the 
% large ROI (~200cm^2). Small ROI values are excluded which have their
% central point outside the large ROI. The percentage
% intensity uniformity is then calculated based on the eligible ROIs with 
% the maximum and minimum means. The results are visualised.

function PIU = ACR_Uniformity(img_ACR,obj_ACR,options)

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_unif = squeeze(double(img_ACR(:,:,7,:))); % if yes, only process the first
else
    img_unif = double(img_ACR(:,:,7));
end

for q = 1:size(img_ACR,4)
    res_ACR = ACR_RetrievePixelSpacing(obj_ACR(q));

    r_img = ceil(80/res_ACR(1)); % equivalent pixel radius
    r_small = ceil(sqrt(100/pi)/res_ACR(1)); % equivalent pixel radius for small ROI
    d_void = ceil(5/res_ACR(1)); % distance from top of phantom to end of void in mm

    % Find centroid
    centroid = ACR_Centroid(img_ACR(:,:,:,q),obj_ACR(q)); % determine centroid

    % ROI
    [img_cols, img_rows] = meshgrid(1:size(img_unif,1),1:size(img_unif,2)); % create grid
    roi_index = (img_rows - centroid(2) - d_void).^2 + (img_cols - centroid(1)).^2 <= r_img.^2; % create large ROI

    % ACR Percentage Integral Uniformity (PIU)
    roi = zeros(size(roi_index)); % pre-allocate small ROI arrays
    mean_array = zeros(size(roi_index)); % pre-allocate mean of small ROIs

    img_unif_masked = roi_index.*img_unif(:,:,q); % mask uniformity slice with large ROI

    base_mask = (img_rows - centroid(2) - d_void).^2 + (img_cols - centroid(1)).^2 <= r_small.^2;
    [rows,cols] = find(base_mask);
    [lrows,lcols] = find(roi_index);

    % tic
    for k = 1:nnz(roi_index)
        centre = [lcols(k),lrows(k)];
        trans_mask = [rows+centre(2)-centroid(2)-d_void,cols+centre(1)-centroid(1)];
        ind = sub2ind(size(img_unif_masked),trans_mask(:,2),trans_mask(:,1));
        roi = img_unif_masked(ind);
        if nnz(roi) < nnz(base_mask)
            mean_array(lcols(k),lrows(k)) = 0;
        else
            mean_array(lcols(k),lrows(k)) = mean(nonzeros(roi),'all');
        end
    end
    % toc

    sig_max = max(mean_array(:)); % find max signal
    sig_min = min(nonzeros(mean_array(:))); % find min signal

    [max_row, max_col] = find(mean_array == sig_max,1); % find centre of max signal ROI
    [min_row, min_col] = find(mean_array == sig_min,1); % find centre of min signal ROI

    PIU(q) = 100*(1 - (sig_max-sig_min)/(sig_max+sig_min)); % Conventional ACR uniformity

    if strcmp(options.SuppressFigures,'no')
        figure
        imshow(img_unif(:,:,q),[])
        hold on
        plot(r_img*cosd(0:1:360)+centroid(1),r_img*sind(0:1:360)+centroid(2)+d_void)
        plot([max_col min_col],[max_row min_row],'r*')
        plot(r_small*cosd(0:1:360)+max_col,r_small*sind(0:1:360)+max_row,'color','y')
        text(max_col, max_row+floor(10./res_ACR(1)),... T
            sprintf('max = %.1f',sig_max),'color','w','fontsize',10) % label with measured mean
        plot(r_small*cosd(0:1:360)+min_col,r_small*sind(0:1:360)+min_row,'color','y')
        text(min_col, min_row+floor(10./res_ACR(1)),...
            sprintf('min = %.1f',sig_min),'color','w','fontsize',10) % label with measured mean
        % xlabel(['PIU = ' num2str(round(PIU,2)) '%'],'fontweight','bold','fontsize',14)
        % title('Uniformity')
    end

    % NEMA Peak Deviation Non-Uniformity (PDNU)
    smooth_filter = [1 2 1; 2 4 2; 1 2 1]; % low-pass filter
    smooth_img_unif = (1/sum(smooth_filter(:))).*conv2(img_unif(:,:,q),smooth_filter,'same'); % apply smoothing filter
    smooth_img_unif_masked = roi_index.*smooth_img_unif; % mask based on large ROI
    max_smooth = max(smooth_img_unif_masked(roi_index)); % find max
    min_smooth = min(smooth_img_unif_masked(roi_index)); % find min

    PDNU = 100*(max_smooth-min_smooth)/(max_smooth+min_smooth); % calculate PDNU

    % NEMA Normalised Absolute Average Deviation Uniformity (NAADU)
    mean_smooth = mean(smooth_img_unif_masked(roi_index)); % find mean
    diff_smooth = abs(smooth_img_unif_masked(roi_index) - mean_smooth); % find absolute deviation
    NAADU = 100*(1 - (1/(length(diff_smooth)*mean_smooth))*sum(diff_smooth(:)));
end

if isrow(PIU)
    PIU = PIU';
end