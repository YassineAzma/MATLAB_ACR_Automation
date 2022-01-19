%% ACR Slice Position
% by Yassine Azma (Oct 2021)
%
% This script uses the ACR series to take line profiles of the slice position 
% bars ramps. The difference between these two line profiles is then used
% to find the bar offset. The results are visualised.

function dL = ACR_SlicePosition(img_ACR,obj_ACR)
close all
slice_num = [1 11];
for n = 1:2
    figure
    if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
        img_insert = squeeze(double(img_ACR(:,:,slice_num(n),1))); % if yes, only process the first
        waitfor(msgbox('4D array detected. Only processing first axial series.'));
    else
        img_insert = double(img_ACR(:,:,slice_num(n)));
    end

    if isempty(obj_ACR.getAttributeByName('PixelSpacing')) % Multi-frame check
        list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');
        res = list.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
    else
        res = obj_ACR.getAttributeByName('PixelSpacing'); % retrieve ACR in-plane resolution
    end
    
    r = 6.5; % radius of slice position insert in mm
    r_dim = round(6/res(1)); % radius in pixels

    % Find centroid
    thresh = bwareaopen(img_insert>0.2*max(img_insert(:)),500); % threshold and remove unconnected pixels
    bhull = bwconvhull(thresh); % create convex hull image
    centroid = ACR_Centroid(img_ACR);

    first_row = find(sum(bhull,1)>0,1,'first')-1;  % Bottom dot in image #1
    last_row = find(sum(bhull,2)>0,1,'first')-1; % Top dot in image #1

    line_prof = [];
    interp_line_prof = [];
    x = [r_dim*[-0.5;0.5]+centroid(1)+1]';

    imshow(img_insert,[],'InitialMagnification',400)
    hold on
    plot([x(1) x(1)],[first_row+round(10/res(1)) centroid(2)-round(50/res(1))],'b')
    plot([x(2) x(2)],[first_row+round(10/res(1)) centroid(2)-round(50/res(1))],'r')
    hold off

    answer = questdlg('Is the current LR position of the line profiles suitable?',...
        'Line Profile Position Confirmation',...
        'Yes','No','No');

    while strcmp(answer,'No')
        switch answer
            case 'No'
                msgbox('Zoom to the bars and press any key.','modal');
                zoom on
                pause()
                zoom off
                msgbox('Using the crosshair, select the correct LR position of the line profiles.','modal');
                [x,~] = ginput(2);

                x = sort(x);
                imshow(img_insert,[],'InitialMagnification',400)
                hold on
                plot([x(1) x(1)],[first_row+round(10/res(1)) centroid(2)-round(50/res(1))],'b')
                plot([x(2) x(2)],[first_row+round(10/res(1)) centroid(2)-round(50/res(1))],'r')

                answer = questdlg('Is the current LR position of the line profiles suitable?',...
                    'Line Profile Position Confirmation',...
                    'Yes','No','No');
            case 'Yes'
                break
        end
    end

    % Draw line profiles across bars
    for k = 1:2
        line_prof(:,k) = improfile(img_insert,[x(k) x(k)],[first_row+round(10/res(1)) centroid(2)-round(50/res(1))]);
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

    [pk, pk_index] = findpeaks(delta_abs,'SortStr','descend','NPeaks',1); % find peak of delta

    inflect_1 = findchangepts(interp_line_prof(:,1),'MaxNumChanges',2,'Statistic','linear'); % find inflection points
    inflect_2 = findchangepts(interp_line_prof(:,2),'MaxNumChanges',2,'Statistic','linear'); % find inflection points

    min_lag = floor(-1.5*abs(mean(inflect_1 - inflect_2))); % minimum shift
    max_lag = ceil(1.5*abs(mean(inflect_1 - inflect_2))); % maximum shift

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

    err = rms(difference,1); % find rms error of each difference profile

    if pos == 1
        shift = -lag(err==min(err(err>0))); % find lag with smallest rms
    else
        shift = lag(err==min(err(err>0)));
    end

    dL(n) = pos*abs(shift)*(1/interp_factor)*res(2); % calculate bar length difference

    subplot(2,2,[1,3])
    imshow(img_insert,[])
    hold on
    plot([x(1) x(1)],[first_row+round(10/res(1)) centroid(2)-round(50/res(1))],'b')
    plot([x(2) x(2)],[first_row+round(10/res(1)) centroid(2)-round(50/res(1))],'r')
    hold off
    % xlabel(['\DeltaL = ' num2str(round(dL,2)) 'mm'],'fontweight','bold','fontsize',14)
    title('Image')

    subplot(2,2,2)
    plot((1/interp_factor)*[1:length(interp_line_prof(:,1))].*res(2),interp_line_prof(:,1),'b')
    hold on
    plot((1/interp_factor)*[1:length(interp_line_prof(:,2))].*res(2),interp_line_prof(:,2),'r')
    xlabel('Distance (mm)')
    ylabel('Intensity')
    grid on
    title('Line Profiles')

    subplot(2,2,4)
    plot((1/interp_factor)*[1:length(interp_line_prof(:,1))].*res(2),interp_line_prof(:,1),'b')
    hold on
    plot((1/interp_factor)*[1:length(interp_line_prof(:,2))].*res(2),circshift(interp_line_prof(:,2),pos*shift),'r')
    xline((1/interp_factor)*left_index.*res(1),'--')
    xline((1/interp_factor)*right_index.*res(1),'--')
    xlabel('Distance (mm)')
    ylabel('Intensity')
    grid on
    title('Shifted')
    % sgtitle('Slice Position')
end