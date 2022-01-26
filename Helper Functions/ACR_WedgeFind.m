function [x,y] = ACR_WedgeFind(img,mask,res_ACR)

%% X LOCATION 
investigate_region = ceil(35/res_ACR(1));

if ~mod(investigate_region,2)
    investigate_region = investigate_region + 1;
end

w_point = find(sum(mask,1)>0,1,'first')-1; % westmost point
e_point = find(sum(mask,1)>0,1,'last')-1; % eastmost point
n_point = find(sum(mask,2)>0,1,'first')-1; % northmost point

for k = 1:investigate_region
    y_loc = n_point+k;
    t(:,k) = mask(y_loc,w_point:e_point)'.*improfile(img,[w_point e_point],y_loc*[1 1]);
end

mean_x_profile = mean(t,2);
abs_diff_x_profile = abs(diff(mean_x_profile));
[~,x_locs] = findpeaks(abs_diff_x_profile,'NPeaks',2,'SortStr','descend');
x_locs = w_point + x_locs-1;

width_pts = [x_locs(1) x_locs(2)];
width = max(width_pts)-min(width_pts);

x = round([min(width_pts)+0.33*width, max(width_pts)-0.33*width]);
%% Y LOCATION
investigate_region = ceil(20/res_ACR(2)); % Width to investigate

end_point = n_point + round(50/res_ACR(2));

if ~mod(investigate_region,2)
    investigate_region = investigate_region + 1;
end

for k = 1:investigate_region
    x_loc = k-floor(investigate_region/2)+floor(mean(x));
    c(:,k) = mask(n_point:end_point,x_loc).*improfile(img,x_loc*[1 1],[n_point end_point]);
end

mean_y_profile = mean(c,2);
abs_diff_y_profile = abs(diff(mean_y_profile));

[~,y_locs] = findpeaks(abs_diff_y_profile,'NPeaks',2,'SortStr','descend');

y_locs = n_point + y_locs-1;

if y_locs(2)-y_locs(1) < 5/res_ACR(2)
    y = n_point + round(10/res_ACR(2));
else
    y = round(min(y_locs)+0.4*abs(diff(y_locs)));
end
%% LENGTH

dist_to_y = abs(n_point - y(1))*res_ACR(2); % in mm
y(2) = round(y(1)+(47 - dist_to_y)/res_ACR(2));

imshow(img,[])
hold on
plot([x(1),x(1)],[y(1),y(2)],'-b')
plot([x(2),x(2)],[y(1),y(2)],'-r')
