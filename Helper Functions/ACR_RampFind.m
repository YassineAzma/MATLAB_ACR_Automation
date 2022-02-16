function [x,y] = ACR_RampFind(img,centroid,res_ACR)

%% X LOCATION 
investigate_region = ceil(5.5/res_ACR(2));

if ~mod(investigate_region,2)
    investigate_region = investigate_region + 1;
end

for k = 1:investigate_region
    y_loc = centroid(2)+k;
    t(:,k) = improfile(img,[1 size(img,2)],y_loc*[1 1]);
end

mean_x_profile = mean(t,2);
abs_diff_x_profile = abs(diff(mean_x_profile));
[~,xlocs] = findpeaks(abs_diff_x_profile,'NPeaks',4,'SortStr','descend');
xlocs = sort(xlocs)-1;

width_pts = [xlocs(2) xlocs(3)];
width = max(width_pts)-min(width_pts);

x = round([min(width_pts)+0.2*width, max(width_pts)-0.2*width]);

%% Y LOCATION

c = improfile(img,[centroid(1) centroid(1)],centroid(2)+[-2*investigate_region 2*investigate_region]);

abs_diff_y_profile = abs(diff(c));
[~,ylocs] = findpeaks(abs_diff_y_profile,'NPeaks',2,'SortStr','descend');

height_pts = centroid(2)-2*investigate_region-1+ylocs;
height = max(height_pts)-min(height_pts);

y = round([max(height_pts)-0.3*height min(height_pts)+0.3*height]);

% imshow(img,[])
% hold on
% plot([x(1),x(2)],[y(1),y(1)],'-b')
% plot([x(1),x(2)],[y(2),y(2)],'-r')