%% ACR Slice Thickness
% by Yassine Azma (Oct 2021)
%
% This script uses the ACR series to take line profiles of the ramps within
% the insert. The full-width half-maximum of the respective profiles is then 
% used to find the ramp lengths. The results are visualised.

function dz = ACR_SliceThickness(img_ACR,obj_ACR)
close all

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_insert = squeeze(double(img_ACR(:,:,1,1))); % if yes, only process the first
    waitfor(msgbox('4D array detected. Only processing first axial series.'));
else
    img_insert = double(img_ACR(:,:,1));
end

if isempty(obj_ACR.getAttributeByName('PixelSpacing')) % Multi-frame check
    list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');
    res = list.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
else
    res = obj_ACR.getAttributeByName('PixelSpacing'); % retrieve ACR in-plane resolution
end

if isempty(obj_ACR.getAttributeByName('SliceThickness')) % Multi-frame check
    list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');
    dz_actual = list.Item_1.PixelMeasuresSequence.Item_1.SliceThickness;
else
    dz_actual = obj_ACR.getAttributeByName('SliceThickness'); % retrieve ACR prescribed slice thickness
end

dims = [180; 8]; % dimensions of insert in cm
dims_img = round(dims.*(1./res)); % dimensions of insert in pixels

% Find centroid
bhull = bwmorph(bwconvhull(img_insert>0.3*max(img_insert(:))),'thin',1); % create convex hull image
centroid = floor(regionprops(bhull,'Centroid').Centroid); % determine centroid from convex hull image

% Find rotation of phantom
% TO DO: LINE PROFILES TO FIND ASSYMMETRY?


% initial line profile coordinates through ramp
x = round([centroid(1) - floor(0.6*dims_img(1)/2), centroid(1) + floor(0.6*dims_img(1)/2),... 
    centroid(1) - floor(0.6*dims_img(1)/2), centroid(1) + floor(0.6*dims_img(1)/2)]); % x points
y = round([centroid(2) + floor(dims_img(2)*0.35), centroid(2) + floor(dims_img(2)*0.35),... 
    centroid(2) - floor(dims_img(2)*0.35), centroid(2) - floor(dims_img(2)*0.35)]); % y points

imshow(img_insert,[0 0.1*max(img_insert(:))],'InitialMagnification',400)
hold on
plot([x(1),x(2)],[y(1),y(2)],'b')
plot([x(3),x(4)],[y(3),y(4)],'r')
hold off

answer = questdlg('Is the current AP position of the line profiles suitable?',...
    'Line Profile Position Confirmation',...
    'Yes','No','No');

while strcmp(answer,'No')
    switch answer
        case 'No'
            msgbox('Zoom to the ramps and press any key.','modal');
            zoom on
            pause()
            zoom off
            msgbox('Using the crosshair, select the correct AP position of the line profiles.','modal');
            [~,y] = ginput(2);

            y = [y(1) y(1) y(2) y(2)];
            imshow(img_insert,[0 0.1*max(img_insert(:))],'InitialMagnification',400)
            hold on
            plot([x(1),x(2)],[y(1) y(2)],'b')
            plot([x(3),x(4)],[y(3) y(4)],'r')

            answer = questdlg('Is the current AP position of the line profiles suitable?',...
                'Line Profile Position Confirmation',...
                'Yes','No','No');
        case 'Yes'
            break
    end
end

interp_factor = 5; % interpolation factor for subpixel measurement
ramp_prof = zeros(max(diff(x))*interp_factor,7,2);
temp = 0;

for m = 1:5 % Try different offsets around centre of phantom
    for k = 1:2
        % Prepare line profile data
        for n = 1:3 % sum line profiles around centre
            line_prof = [];
            line_prof = improfile(img_insert,[x(2*k-1) x(2*k)],-3+m+(n-2)+[y(2*k-1) y(2*k)]);
            line_prof = interp1(line_prof,0:(1/interp_factor):length(line_prof),'pchip');
            baseline = min(line_prof); % baseline
            line_prof = line_prof - baseline; % remove baseline
            temp = temp+line_prof;
        end

        line_prof = temp/n;
        % Find FWHM of line profile
        [pk(m,k), pk_index(m,k)] = maxk(line_prof,1);
        [~, left_index(m,k)] = min(abs(line_prof(1:pk_index(m,k)) - 0.5*pk(m,k)));
        [~, right_index(m,k)] = min(abs(line_prof(pk_index(m,k):end) - 0.5*pk(m,k)));
        right_index(m,k) = right_index(m,k) + pk_index(m,k);

        % Calculate ramp length
        index(:,m,k) = [left_index(m,k) right_index(m,k)];
        ramp_length(m,k) = (1/interp_factor)*(diff(index(:,m,k))).*res(1);
        ramp_prof(1:length(line_prof),m,k) = line_prof;

        temp = 0;
    end
end

dz_list = 0.2*(ramp_length(:,1).*ramp_length(:,2)./(ramp_length(:,1)+ramp_length(:,2))); % calculate slice thicknesses

[~,closestIndex] = min(abs(dz_list-dz_actual)); % find closest to desired 
dz = dz_list(closestIndex); % select slice thickness closest to desired

% figure
subplot(2,2,[1,3])
imshow(img_insert,[])
hold on
plot(-2+closestIndex+[x(1) x(2)],[y(1) y(2)],'b')
plot(-2+closestIndex+[x(3) x(4)],[y(3) y(4)],'r')
hold off
% xlabel(['\Deltaz = ' num2str(round(dz,2)) 'mm'],'fontweight','bold','fontsize',14)
title('Image')

subplot(2,2,2)
plot([1:length(ramp_prof(:,closestIndex,2))].*res(1),ramp_prof(:,closestIndex,2),'r')
hold on
plot(res(2).*[pk_index(closestIndex,2) pk_index(closestIndex,2)],[0 pk(closestIndex,2)],'--rx','MarkerIndices',2)
plot(res(2).*index(:,closestIndex,2),0.5*[pk(closestIndex,2) pk(closestIndex,2)],'--r*')
text(res(2).*pk_index(closestIndex,2), 0.6*pk(closestIndex,2),...
    sprintf('L = %.1fmm',ramp_length(closestIndex,2)),'color','k') % label with measured distance
xlabel('Distance (mm)')
ylabel('Intensity')
grid on
title('Upper Ramp')

subplot(2,2,4)
plot([1:length(ramp_prof(:,closestIndex,1))].*res(1),ramp_prof(:,closestIndex,1),'b')
hold on
plot(res(1).*[pk_index(closestIndex,1) pk_index(closestIndex,1)],[0 pk(closestIndex,1)],'--bx','MarkerIndices',2)
plot(res(1).*index(:,closestIndex,1),0.5*[pk(closestIndex,1) pk(closestIndex,1)],'--b*')
text(res(1).*pk_index(closestIndex,1), 0.6*pk(closestIndex,1),...
    sprintf('L = %.1fmm',ramp_length(closestIndex,1)),'color','k') % label with measured distance
xlabel('Distance (mm)')
ylabel('Intensity')
grid on
title('Lower Ramp')
% sgtitle('Slice Thickness')