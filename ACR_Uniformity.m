%% ACR Uniformity
% by Yassine Azma (Oct 2021)
%
% This script uses the ACR series to place small ROIs (~1cm^2) within the 
% large ROI (~200cm^2). Small ROI values are excluded which have their
% central point outside the large ROI. The percentage
% intensity uniformity is then calculated based on the eligible ROIs with 
% the maximum and minimum means. The results are visualised.

function PIU = ACR_Uniformity(img_ACR,obj_ACR)
close all

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_unif = double(img_ACR(:,:,7,1)); % if yes, only process the first
    waitfor(msgbox('4D array detected. Only processing first axial series.'));
else
    img_unif = double(img_ACR(:,:,7));
end

if isempty(obj_ACR.getAttributeByName('PixelSpacing')) % Multi-frame check
    list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');
    res = list.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
else
    res = obj_ACR.getAttributeByName('PixelSpacing'); % retrieve ACR in-plane resolution
end

r = 8; % for ~201cm^2 circular ROI
r_img = ceil(10*r./res(1)); % equivalent pixel radius
r_small = sqrt(100/pi)/res(1); % equivalent pixel radius for small ROI
d_void = 5; % distance from top of phantom to end of void in mm
thresh = 40;

% Find centroid
bhull = bwconvhull(img_unif>thresh/100*max(img_unif(:))); % create convex hull image
centroid = floor(regionprops(bhull,'Centroid').Centroid); % determine centroid from convex hull image

imshow(img_unif,[],'InitialMagnification',400)
hold on
plot(centroid(1),centroid(2),'rx','MarkerSize',12)

answer = questdlg('Is the current centroid location correct?',...
    'Centroid Location Confirmation',...
    'Yes','No','No');

while strcmp(answer,'No')
    switch answer
        case 'No'
            prompt = {['Enter new threshold [0-100] (current = ' num2str(thresh) '):']};
            dlgtitle = 'Input';
            dims = [1 50];
            definput = {'20'};
            thresh = inputdlg(prompt,dlgtitle,dims,definput);
            thresh = str2double(thresh);

            bhull = bwconvhull(img_unif>thresh/100*max(img_unif(:))); % create binary image
            centroid = floor(regionprops(bhull,'Centroid').Centroid); % determine centroid from convex hull image

            imshow(img_unif,[],'InitialMagnification',400)
            hold on
            plot(centroid(1),centroid(2),'rx','MarkerSize',12)

            pause(1)
            answer = questdlg('Is the current centroid location suitable?',...
                'Centroid Location Confirmation',...
                'Yes','No','No');
        case 'Yes'
            break
    end
end

% ROI
[img_cols, img_rows] = meshgrid(1:size(img_unif,1),1:size(img_unif,2)); % create grid 
roi_index = (img_rows - centroid(2) - d_void/res(1)).^2 + (img_cols - centroid(1)).^2 <= r_img.^2; % create large ROI 

% ACR Percentage Integral Uniformity (PIU)
roi = zeros(size(roi_index)); % pre-allocate small ROI arrays
mean_val = zeros(size(roi_index)); % pre-allocate mean of small ROIs

img_unif_masked = roi_index.*img_unif; % mask uniformity slice with large ROI

f = waitbar(0,'1','Name','Placing Small ROIs...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for i = 1:size(img_unif,1)
    waitbar(i/size(img_unif,1),f,sprintf('Simulating Column %.0f of %.0f',i,size(img_unif,1)))
    for j = 1:size(img_unif,2)
        if roi_index(i,j) == 0
            mean_val(i,j) = 0; % don't create small ROI if centre outside large ROI
        else     
            roi = (img_rows - i).^2 + (img_cols - j).^2 <= r_small.^2; % create ~1cm^2 small circular ROI
            roi_vals = img_unif_masked(roi); % retrieve small ROI values
            if length(nonzeros(roi_vals)) < length(find(roi))
                mean_val(i,j) = 0; % don't create small ROI if centre outside large ROI
            else
                mean_val(i,j) = mean(nonzeros(roi_vals)); % take mean
            end
        end
    end
end
delete(f)

sig_max = max(nonzeros(mean_val(:))); % find max signal
sig_min = min(nonzeros(mean_val(:))); % find min signal

[max_row, max_col] = find(mean_val == sig_max,1); % find centre of max signal ROI
[min_row, min_col] = find(mean_val == sig_min,1); % find centre of min signal ROI

PIU = 100*(1 - (sig_max-sig_min)/(sig_max+sig_min)); % Conventional ACR uniformity

% NEMA Peak Deviation Non-Uniformity (PDNU)
smooth_filter = [1 2 1; 2 4 2; 1 2 1]; % low-pass filter
smooth_img_unif = (1/sum(smooth_filter(:))).*conv2(img_unif,smooth_filter,'same'); % apply smoothing filter
smooth_img_unif_masked = roi_index.*smooth_img_unif; % mask based on large ROI
max_smooth = max(smooth_img_unif_masked(roi_index)); % find max
min_smooth = min(smooth_img_unif_masked(roi_index)); % find min

PDNU = 100*(max_smooth-min_smooth)/(max_smooth+min_smooth); % calculate PDNU

% NEMA Normalised Absolute Average Deviation Uniformity (NAADU)
mean_smooth = mean(smooth_img_unif_masked(roi_index)); % find mean
diff_smooth = abs(smooth_img_unif_masked(roi_index) - mean_smooth); % find absolute deviation
NAADU = 100*(1 - (1/(length(diff_smooth)*mean_smooth))*sum(diff_smooth(:))); 

imshow(img_unif,[])
hold on
plot(r_img*cosd(0:1:360)+centroid(1),r_img*sind(0:1:360)+centroid(2)+d_void/res(1))
plot([max_col min_col],[max_row min_row],'r*')
plot(r_small*cosd(0:1:360)+max_col,r_small*sind(0:1:360)+max_row,'color','y')
text(max_col, max_row+floor(10./res(1)),...
    sprintf('max = %.1f',sig_max),'color','w','fontsize',10) % label with measured mean
plot(r_small*cosd(0:1:360)+min_col,r_small*sind(0:1:360)+min_row,'color','y')
text(min_col, min_row+floor(10./res(1)),...
    sprintf('min = %.1f',sig_min),'color','w','fontsize',10) % label with measured mean
% xlabel(['PIU = ' num2str(round(PIU,2)) '%'],'fontweight','bold','fontsize',14)
% title('Uniformity')