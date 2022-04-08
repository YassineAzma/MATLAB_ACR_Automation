%% ACR SNR
% by Yassine Azma (Oct 2021)
%
% This script uses the ACR series to calculate the SNR. 
% The Rician noise correction can be included in the options structure. 
% The results are visualised.

function SNR = ACR_SNR(img_ACR,obj_ACR,options)

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_SNR = double(img_ACR(:,:,7,:)); % if yes, only process the first
else
    img_SNR = double(img_ACR(:,:,7));
end

for k = 1:size(img_ACR,4)
    res_ACR = ACR_RetrievePixelSpacing(obj_ACR(k));

    r = 8; % for ~201cm^2 circular ROI
    r_img = ceil(10*r./res_ACR(1)); % equivalent pixel radius

    % Find centroid
    centroid = ACR_Centroid(img_ACR(:,:,:,k),obj_ACR(k)); % determine centroid
    temp = img_SNR(:,:,k);
    mask = ACR_Threshold(temp,res_ACR,centroid);

    % ROI
    [img_cols, img_rows] = meshgrid(1:size(temp,1),1:size(temp,2)); % create grid
    roi_index = (img_rows - centroid(2) - 5).^2 + (img_cols - centroid(1)).^2 <= r_img.^2; % create large ROI

    w_point = 0.75*find(sum(mask,1)>0,1,'first')-1; % westmost point
    e_point = 0.5*w_point + find(sum(mask,1)>0,1,'last')-1; % eastmost point
    n_point = 0.75*find(sum(mask,2)>0,1,'first')-1; % northmost point
    s_point = 0.5*n_point + find(sum(mask,2)>0,1,'last')-1; % southmost point

    % Noise ROIs
    nw_noise_index = (img_rows - w_point).^2 + (img_cols - n_point).^2 <= ceil(10./res_ACR(1)).^2; % create noise ROI
    ne_noise_index = (img_rows - n_point).^2 + (img_cols - e_point).^2 <= ceil(10./res_ACR(1)).^2; % create noise ROI
    se_noise_index = (img_rows - e_point).^2 + (img_cols - s_point).^2 <= ceil(10./res_ACR(1)).^2; % create noise ROI
    sw_noise_index = (img_rows - s_point).^2 + (img_cols - w_point).^2 <= ceil(10./res_ACR(1)).^2; % create noise ROI

    sig_mean = mean(nonzeros(temp(roi_index))); % mean signal in large ROI
    nw_noise_std = std(nonzeros(temp(nw_noise_index))); % std of NW noise ROI
    ne_noise_std = std(nonzeros(temp(ne_noise_index))); % std of NE noise ROI
    se_noise_std = std(nonzeros(temp(se_noise_index))); % std of SE noise ROI
    sw_noise_std = std(nonzeros(temp(sw_noise_index))); % std of SW noise ROI

    SNR(k) = sig_mean/mean([nw_noise_std,ne_noise_std,se_noise_std,sw_noise_std]);

    if strcmp(options.RicianNoiseCorrection,'yes')
        SNR(k) = sqrt((4-pi)/2)*SNR(k);
    end

    if strcmp(options.SuppressFigures,'no')
        figure
        imshow(img_SNR,[0 0.1*max(temp(:))])
        hold on
        plot(r_img*cosd(0:1:360)+centroid(1),r_img*sind(0:1:360)+centroid(2)+5)
        text(centroid(1)-3*floor(10./res_ACR(1)), centroid(2)+floor(10./res_ACR(1)),...
            sprintf('mean = %.1f',sig_mean),'color','w','fontsize',8) % label with measured mean

        plot(ceil(10./res_ACR(1))*cosd(0:1:360)+w_point,ceil(10./res_ACR(1))*sind(0:1:360)+n_point,'color','r')
        text(w_point-floor(10./res_ACR(1)), n_point+floor(10./res_ACR(1)),...
            sprintf('std = %.1f',nw_noise_std),'color','w','fontsize',8) % label with measured std

        plot(ceil(10./res_ACR(1))*cosd(0:1:360)+e_point,ceil(10./res_ACR(1))*sind(0:1:360)+n_point,'color','r')
        text(e_point-floor(10./res_ACR(1)), n_point+floor(10./res_ACR(1)),...
            sprintf('std = %.1f',ne_noise_std),'color','w','fontsize',8) % label with measured std

        plot(ceil(10./res_ACR(1))*cosd(0:1:360)+e_point,ceil(10./res_ACR(1))*sind(0:1:360)+s_point,'color','r')
        text(e_point-floor(10./res_ACR(1)), s_point+floor(10./res_ACR(1)),...
            sprintf('std = %.1f',se_noise_std),'color','w','fontsize',8) % label with measured std

        plot(ceil(10./res_ACR(1))*cosd(0:1:360)+w_point,ceil(10./res_ACR(1))*sind(0:1:360)+s_point,'color','r')
        text(w_point-floor(10./res_ACR(1)), s_point+floor(10./res_ACR(1)),...
            sprintf('std = %.1f',sw_noise_std),'color','w','fontsize',8) % label with measured std

        % xlabel(['SNR = ' num2str(round(SNR,1))],'fontweight','bold','fontsize',14)
        % title('Signal to Noise Ratio')
    end
end

%% Subtraction SNR

if size(img_ACR,4) > 1
    temp = img_SNR(:,:,1);
    sig_mean_1 = mean(nonzeros(temp(roi_index))); % mean signal in large ROI
    temp = img_SNR(:,:,2);
    sig_mean_2 = mean(nonzeros(temp(roi_index)));

    sig_mean = 0.5*(sig_mean_1+sig_mean_2);
    noise_std = (1/sqrt(2))*std(nonzeros(img_SNR_sub(roi_index))); % std of noise ROI

    SNR = sig_mean/noise_std;

    if strcmp(options.SuppressFigures,'no')
        figure
        subplot(1,2,1)
        imshow(img_SNR(:,:,1),[])
        hold on
        plot(r_img*cosd(0:1:360)+centroid(1),r_img*sind(0:1:360)+centroid(2)+5)
        text(centroid(1)-3*floor(10./res_ACR(1)), centroid(2)+floor(10./res_ACR(1)),...
            sprintf('mean = %.1f',sig_mean),'color','w','fontsize',8) % label with averaged mean

        subplot(1,2,2)
        imshow(img_SNR_sub,[])
        hold on
        plot(r_img*cosd(0:1:360)+centroid(1),r_img*sind(0:1:360)+centroid(2)+5)
        text(centroid(1)-3*floor(10./res_ACR(1)), centroid(2)+floor(10./res_ACR(1)),...
            sprintf('std = %.1f',noise_std),'color','w','fontsize',8) % label with measured std

        % xlabel(['SNR = ' num2str(round(SNR,1))],'fontweight','bold','fontsize',14)
        % title('Signal to Noise Ratio')
    end
end
