%% ACR Slice Thickness
% by Yassine Azma (Oct 2021)
%
% This script uses the ACR series to take line profiles of the ramps within
% the insert. The full-width half-maximum of the respective profiles is then 
% used to find the ramp lengths. The results are visualised.

function dz = ACR_SliceThickness(img_ACR,obj_ACR)
close all

% Find rotation of phantom and correct
rot_angle = ACR_FindRotation(img_ACR,obj_ACR);
img_ACR = ACR_RotateImages(img_ACR,rot_angle);

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_insert = squeeze(double(img_ACR(:,:,1,1))); % if yes, only process the first
    waitfor(msgbox('4D array detected. Only processing first axial series.'));
else
    img_insert = double(img_ACR(:,:,1));
end

res_ACR = ACR_RetrievePixelSpacing(obj_ACR);
dz_actual = ACR_RetrieveSliceThickness(obj_ACR);

dims = [180; 8]; % dimensions of insert in cm
dims_img = round(dims.*(1./res_ACR)); % dimensions of insert in pixels

% Find centroid
centroid = ACR_Centroid(img_ACR,obj_ACR); % determine centroid

% initial line profile coordinates through ramp
[x,y] = ACR_RampFind(img_insert,centroid,res_ACR);

x = repmat(x,1,2);
y = repelem(y,2);
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
        ramp_length(m,k) = (1/interp_factor)*(diff(index(:,m,k))).*res_ACR(1);
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
plot([1:length(ramp_prof(:,closestIndex,2))].*res_ACR(1),ramp_prof(:,closestIndex,2),'r')
hold on
plot(res_ACR(2).*[pk_index(closestIndex,2) pk_index(closestIndex,2)],[0 pk(closestIndex,2)],'--rx','MarkerIndices',2)
plot(res_ACR(2).*index(:,closestIndex,2),0.5*[pk(closestIndex,2) pk(closestIndex,2)],'--r*')
text(res_ACR(2).*pk_index(closestIndex,2), 0.6*pk(closestIndex,2),...
    sprintf('L = %.1fmm',ramp_length(closestIndex,2)),'color','k') % label with measured distance
xlabel('Distance (mm)')
ylabel('Intensity')
grid on
title('Upper Ramp')

subplot(2,2,4)
plot([1:length(ramp_prof(:,closestIndex,1))].*res_ACR(1),ramp_prof(:,closestIndex,1),'b')
hold on
plot(res_ACR(1).*[pk_index(closestIndex,1) pk_index(closestIndex,1)],[0 pk(closestIndex,1)],'--bx','MarkerIndices',2)
plot(res_ACR(1).*index(:,closestIndex,1),0.5*[pk(closestIndex,1) pk(closestIndex,1)],'--b*')
text(res_ACR(1).*pk_index(closestIndex,1), 0.6*pk(closestIndex,1),...
    sprintf('L = %.1fmm',ramp_length(closestIndex,1)),'color','k') % label with measured distance
xlabel('Distance (mm)')
ylabel('Intensity')
grid on
title('Lower Ramp')
% sgtitle('Slice Thickness')