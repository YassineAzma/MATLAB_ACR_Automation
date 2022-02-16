%% ACR Slice Position
% by Yassine Azma (Oct 2021)
%
% This script uses the ACR series to take line profiles of the slice position 
% bars ramps. The difference between these two line profiles is then used
% to find the bar offset. The results are visualised.

function dL = ACR_SlicePosition(img_ACR,obj_ACR,options)
close all
slice_num = [1 11];

res_ACR = ACR_RetrievePixelSpacing(obj_ACR);

rot_angle = ACR_FindRotation(img_ACR,obj_ACR);
img_ACR = ACR_RotateImages(img_ACR,rot_angle);
centroid = ACR_Centroid(img_ACR,obj_ACR);
mask = ACR_Threshold(img_ACR,res_ACR,centroid);

if strcmp(options.SliceDisplacementDisplay,'yes')
    corr = 0.5;
else
    corr = 1;
end

for n = 1:2
    figure
    if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
        img_insert = squeeze(double(img_ACR(:,:,slice_num(n),1))); % if yes, only process the first
        waitfor(msgbox('4D array detected. Only processing first axial series.'));
    else
        img_insert = double(img_ACR(:,:,slice_num(n)));
    end

    line_prof = [];
    interp_line_prof = [];
    [x,y] = ACR_WedgeFind(img_insert,mask,res_ACR);

    % Draw line profiles across bars
    for k = 1:2
        line_prof(:,k) = improfile(img_insert,[x(k) x(k)],[y(1) y(2)]);
    end

    interp_factor = 5; % interpolation factor for subpixel measurement
    interp_line_prof = interp1(line_prof,0:(1/interp_factor):length(line_prof),'pchip'); %interpolate line profiles
    interp_line_prof = rmmissing(interp_line_prof);

    delta = diff(interp_line_prof,1,2); % subtract line 1 from line 2

    delta_abs = abs(delta); % take absolute

    if max(-delta) > max(delta) % if true, the right bar must be longer than the left, dL is positive as defined by the guidance
        pos = 1;
    else
        pos = -1;
    end

    inflect_1 = findchangepts(interp_line_prof(:,1),'MaxNumChanges',2,'Statistic','linear'); % find inflection points
    inflect_2 = findchangepts(interp_line_prof(:,2),'MaxNumChanges',2,'Statistic','linear'); % find inflection points

    min_lag = floor(-2*abs(mean(inflect_1 - inflect_2))); % minimum shift
    max_lag = ceil(2*abs(mean(inflect_1 - inflect_2))); % maximum shift

    left_index = min([inflect_1(:);inflect_2(:)],[],'all'); % leftmost inflection point
    right_index = max([inflect_1(:);inflect_2(:)],[],'all'); % rightmost inflection point

    difference = [];
    lag = [];

    for k = min_lag:max_lag
        ind = k-min_lag+1;
        difference(:,ind) = interp_line_prof(left_index:right_index,2)...
            - circshift(interp_line_prof(left_index:right_index,1),k); % subtract lagging/leading line profile 1 from line profile 2
        discontinuity(ind) = findchangepts(difference(:,ind));
        if k < 0
            difference(discontinuity(ind):end,ind) = 0; % set all values before/after the discontinuity (dependent on lag sign) to 0
        else
            difference(1:discontinuity(ind),ind) = 0;
        end
        lag(ind) = k;
    end

%     err = rms(difference,1); % find rms error of each difference profile
    err = mean(difference,1);

    ignore_ind = find(lag==0);
    err_ignore = err(1:end ~= ignore_ind);
    if pos == 1
        temp = err(err_ignore==min(err_ignore(err_ignore>0)));
        shift = -lag(err==temp); % find lag with smallest rms
    else
        temp = err(err_ignore==min(err_ignore(err_ignore>0)));
        shift = lag(err==temp);
    end

    dL(n) = corr*pos*abs(shift)*(1/interp_factor)*res_ACR(2); % calculate bar length difference

    subplot(2,2,[1,3])
    imshow(img_insert,[])
    hold on
    plot([x(1) x(1)],[y(1) y(2)],'b')
    plot([x(2) x(2)],[y(1) y(2)],'r')
    hold off
    % xlabel(['\DeltaL = ' num2str(round(dL,2)) 'mm'],'fontweight','bold','fontsize',14)
    title('Image')

    subplot(2,2,2)
    plot((1/interp_factor)*[1:length(interp_line_prof(:,1))].*res_ACR(2),interp_line_prof(:,1),'b')
    hold on
    plot((1/interp_factor)*[1:length(interp_line_prof(:,2))].*res_ACR(2),interp_line_prof(:,2),'r')
    xlabel('Distance (mm)')
    ylabel('Intensity')
    grid on
    title('Line Profiles')

    subplot(2,2,4)
    plot((1/interp_factor)*[1:length(interp_line_prof(:,1))].*res_ACR(2),interp_line_prof(:,1),'b')
    hold on
    plot((1/interp_factor)*[1:length(interp_line_prof(:,2))].*res_ACR(2),circshift(interp_line_prof(:,2),pos*shift),'r')
    xline((1/interp_factor)*left_index.*res_ACR(1),'--')
    xline((1/interp_factor)*right_index.*res_ACR(1),'--')
    xlabel('Distance (mm)')
    ylabel('Intensity')
    grid on
    title('Shifted')
    % sgtitle('Slice Position')
end