%% ACR High Contrast Spatial Resolution
% by Yassine Azma (Nov 2021)
% 
% This script ...

function resolvable = ACR_SpatialResolution(img_ACR,obj_ACR)
%% Masking
close all

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_insert = squeeze(double(img_ACR(:,:,1,1))); % if yes, only process the first
    waitfor(msgbox('4D array detected. Only processing first axial series.'));
else
    img_insert = double(img_ACR(:,:,1));
end

if isempty(obj_ACR.getAttributeByName('PixelSpacing')) % Multi-frame check
    list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');
    res_ACR = list.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing; % retrieve ACR in-plane resolution
else
    res_ACR = obj_ACR.getAttributeByName('PixelSpacing'); % retrieve ACR in-plane resolution
end

thresh = bwareaopen(img_insert > 0.15*max(img_insert(:)),5); % threshold and remove small unconnected pixels
centroid_mask = bwconvhull(thresh); % convex hull image
centroid = ACR_Centroid(img_ACR); % find centroid

mask = img_insert < 0.2*max(img_insert(:)); % threshold
cropped_mask = mask(round(centroid(2)):end,:); % crop image based on centroid
cropped_mask = bwareaopen(imfill(cropped_mask,'holes'),50); % fill holes and remove unconnected pixels
cropped_mask = imclearborder(cropped_mask,18); % clear border

w_point = find(sum(cropped_mask,1)>0,1,'first')-1 + round(32/res_ACR(1)); % westmost point
e_point = find(sum(cropped_mask,1)>0,1,'last')-1 - round(10/res_ACR(1)); % eastmost point
c = improfile(cropped_mask,[w_point+1 w_point+1],[1 size(cropped_mask,1)]);
if isempty(find(c,1))
    c = improfile(cropped_mask,[w_point-1 w_point-1],[1 size(cropped_mask,1)]);
end

extent = find(c); % find non-zeros 
n_point = extent(1) + round(8.6/res_ACR(2)); % northmost point of resolution insert
s_point = extent(end) - round(8.6/res_ACR(2)); % southmost point of resolution insert

hole_size_ul = [2.2/(res_ACR(1)) 2/(res_ACR(1)) 1.8/(res_ACR(1))]; % size in pixels
hole_size_lr = [2.2/(res_ACR(2)) 2/(res_ACR(2)) 1.8/(res_ACR(2))]; % size in pixels
 
%% 1.1 Horizontal UL Array
line_prof = [];
for m = 1:(s_point-n_point)/2
    line_prof(:,m) = improfile(img_ACR(round(centroid(2)):end,:,1),[w_point w_point + round(15/res_ACR(1))],[n_point+m n_point+m]);
    baseline(m) = min(line_prof(:,m));
    temp_pks = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    if length(temp_pks) < 4
        pks(m,:) = zeros(1,4);
        locs(m,:) = zeros(1,4);
        w(m,:) = zeros(1,4);
        prom(m,:) = zeros(1,4);
    else
        [pks(m,:),locs(m,:),w(m,:),prom(m,:)] = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    end
end

pk_spacing = sort(locs,2);
diff_pk_spacing = diff(pk_spacing,1,2);
spacing_test = diff_pk_spacing <= ceil(hole_size_ul(1)); % peaks should be relatively close
spacing_pass = sum(spacing_test,2) >= 1; % pass if at least two peaks have an adequately spaced neighbour
w_test = 0.25*hole_size_ul(1) < w & w < 1.75*hole_size_ul(1); % peaks should have comparable width to feature
pk_test = pks > 0.5*max(line_prof(:)); % pass if peaks satisfy threshold
prom_test = prom > 0.2*max(line_prof(:)); % pass if prominences satisfy threshold

ul_resolv_11 = find(sum(spacing_pass.*w_test.*pk_test.*prom_test,2) >= 3); % output line profiles that pass
temp = ul_resolv_11; % set temp array

if ~isempty(find(diff(temp)>1))
    grouping = find(diff(temp)>1); % find clusters of line profiles that pass

    group = temp(1:grouping(1)); % find first cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    ul_resolv_11(1) = group(ind(1)); % take greatest peak from the cluster
    if length(grouping) >= 1
        for m = 1:length(grouping)-1
            group = temp(grouping(m)+1:grouping(m+1)); % find second and third clusters
            [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
            ul_resolv_11(m+1) = group(ind(1)); % take greatest peak from each cluster
        end
    end
    group = temp(grouping(end)+1:end); % find last cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    ul_resolv_11(length(grouping)+1) = group(ind(1)); % take greatest peak from the cluster

    ul_resolv_11(length(grouping)+2:end) = [];
    ul_resolv_11 = unique(ul_resolv_11);

    if length(ul_resolv_11) < 4 % check if large cluster has been overlooked
        unused = temp(~ismember(temp,ul_resolv_11)); % find passing line profiles not used
        cluster_check = sum(abs(unused - ul_resolv_11') > 1,2) == length(ul_resolv_11); % check which do not neighbour a used line profile
        [~,ind] = sort(sum(pks(unused(cluster_check),:),2),'descend'); % sort by peak height
        if ~isempty(ind)
            ul_resolv_11(end+1) = unused(ind(1)); % take greatest peak
            ul_resolv_11 = sort(ul_resolv_11);
        end
    end
elseif ~isempty(ul_resolv_11) && length(ul_resolv_11) > 4
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    ul_resolv_11 = temp(ind(1:4)); % take greatest peak
    ul_resolv_11 = sort(ul_resolv_11);
elseif ~isempty(temp) 
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    ul_resolv_11 = temp(ind(1)); % take greatest peak
end

if isrow(ul_resolv_11)
    ul_resolv_11 = ul_resolv_11';
end
%% 1.1 Vertical LR Array
line_prof = [];
for m = 1:(s_point - n_point)/2
    line_prof(:,m) = improfile(img_ACR(round(centroid(2)):end,:,1),[w_point+round(10/res_ACR(2))+m w_point+round(10/res_ACR(2))+m],[n_point s_point]);
    baseline(m) = min(line_prof(:,m));
    temp_pks = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    if length(temp_pks) < 4
        pks(m,:) = zeros(1,4);
        locs(m,:) = zeros(1,4);
        w(m,:) = zeros(1,4);
        prom(m,:) = zeros(1,4);
    else
        [pks(m,:),locs(m,:),w(m,:),prom(m,:)] = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    end
end

pk_spacing = sort(locs,2);
diff_pk_spacing = diff(pk_spacing,1,2);
spacing_test = diff_pk_spacing <= ceil(hole_size_lr(1)); % peaks should be relatively close
spacing_pass = sum(spacing_test,2) >= 1; % pass if at least two peaks have an adequately spaced neighbour
w_test = 0.25*hole_size_lr(1) < w & w < 1.75*hole_size_lr(1); % peaks should have reasonable FWHM compared to feature width
pk_test = pks > 0.5*max(line_prof(:)); % pass if peaks satisfy threshold
prom_test = prom > 0.2*max(line_prof(:)); % pass if prominences satisfy threshold

lr_resolv_11 = find(sum(spacing_pass.*w_test.*pk_test.*prom_test,2) >= 3); % output line profiles that pass
temp = lr_resolv_11; % set temp array

if ~isempty(find(diff(temp)>1))
    grouping = find(diff(temp)>1); % find clusters of line profiles that pass

    group = temp(1:grouping(1)); % find first cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    lr_resolv_11(1) = group(ind(1)); % take greatest peak from the cluster
    if length(grouping) >= 2
        for m = 1:length(grouping)-1
            group = temp(grouping(m)+1:grouping(m+1)); % find second and third clusters
            [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
            lr_resolv_11(m+1) = group(ind(1)); % take greatest peak from each cluster
        end
    end
    group = temp(grouping(end)+1:end); % find last cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    lr_resolv_11(length(grouping)+1) = group(ind(1)); % take greatest peak from the cluster

    lr_resolv_11(length(grouping)+2:end) = [];
    lr_resolv_11 = unique(lr_resolv_11);

    if length(lr_resolv_11) < 4 % check if large cluster has been overlooked
        unused = temp(~ismember(temp,lr_resolv_11)); % find passing line profiles not used
        cluster_check = sum(abs(unused - lr_resolv_11') > 1,2) == length(lr_resolv_11); % check which do not neighbour a used line profile
        [~,ind] = sort(sum(pks(unused(cluster_check),:),2),'descend'); % sort by peak height
        if ~isempty(ind)
            lr_resolv_11(end+1) = unused(ind(1)); % take greatest peak
            lr_resolv_11 = sort(lr_resolv_11);
        end
    end
elseif ~isempty(lr_resolv_11) && length(lr_resolv_11) > 4
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    lr_resolv_11 = temp(ind(1:4)); % take greatest peak
    lr_resolv_11 = sort(lr_resolv_11);
elseif ~isempty(temp) 
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    lr_resolv_11 = temp(ind(1)); % take greatest peak
end

if isrow(lr_resolv_11)
    lr_resolv_11 = lr_resolv_11';
end
%% 1.0 Horizontal UL Array

line_prof = [];
for m = 1:(s_point - n_point)/2
    line_prof(:,m) = improfile(img_ACR(round(centroid(2)):end,:,1),[w_point+round(25/res_ACR(1)) w_point + round(45/res_ACR(1))],[n_point+m n_point+m]);
    baseline(m) = min(line_prof(:,m));
    temp_pks = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    if length(temp_pks) < 4
        pks(m,:) = zeros(1,4);
        locs(m,:) = zeros(1,4);
        w(m,:) = zeros(1,4);
        prom(m,:) = zeros(1,4);
    else
        [pks(m,:),locs(m,:),w(m,:),prom(m,:)] = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    end
end

pk_spacing = sort(locs,2);
diff_pk_spacing = diff(pk_spacing,1,2);
spacing_test = diff_pk_spacing <= ceil(hole_size_ul(2)); % peaks should be relatively close
spacing_pass = sum(spacing_test,2) >= 1; % pass if at least two peaks have an adequately spaced neighbour
w_test = 0.25*hole_size_ul(2) < w & w < 1.75*hole_size_ul(2); % peaks should have comparable width to feature
pk_test = pks > 0.4*max(line_prof(:)); % pass if peaks satisfy threshold
prom_test = prom > 0.2*max(line_prof(:)); % pass if prominences satisfy threshold

ul_resolv_10 = find(sum(spacing_pass.*w_test.*pk_test.*prom_test,2) >= 3); % output line profiles that pass
temp = ul_resolv_10; % set temp array

if ~isempty(find(diff(temp)>1))
    grouping = find(diff(temp)>1); % find clusters of line profiles that pass

    group = temp(1:grouping(1)); % find first cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    ul_resolv_10(1) = group(ind(1)); % take greatest peak from the cluster
    if length(grouping) >= 2
        for m = 1:length(grouping)-1
            group = temp(grouping(m)+1:grouping(m+1)); % find second and third clusters
            [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
            ul_resolv_10(m+1) = group(ind(1)); % take greatest peak from each cluster
        end
    end
    group = temp(grouping(end)+1:end); % find last cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    ul_resolv_10(length(grouping)+1) = group(ind(1)); % take greatest peak from the cluster

    ul_resolv_10(length(grouping)+2:end) = [];
    ul_resolv_10 = unique(ul_resolv_10);

    if length(ul_resolv_10) < 4 % check if large cluster has been overlooked
        unused = temp(~ismember(temp,ul_resolv_10)); % find passing line profiles not used
        cluster_check = sum(abs(unused - ul_resolv_10') > 1,2) == length(ul_resolv_10); % check which do not neighbour a used line profile
        [~,ind] = sort(sum(pks(unused(cluster_check),:),2),'descend'); % sort by peak height
        if ~isempty(ind)
            ul_resolv_10(end+1) = unused(ind(1)); % take greatest peak
            ul_resolv_10 = sort(ul_resolv_10);
        end
    end
elseif ~isempty(ul_resolv_10) && length(ul_resolv_10) > 4
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    ul_resolv_10 = temp(ind(1:4)); % take greatest peak
    ul_resolv_10 = sort(ul_resolv_10);
elseif ~isempty(temp) 
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    ul_resolv_10 = temp(ind(1)); % take greatest peak
end

if isrow(ul_resolv_10)
    ul_resolv_10 = ul_resolv_10';
end
%% 1.0 Vertical LR Array

line_prof = [];
for m = 1:(s_point - n_point)/2
    line_prof(:,m) = improfile(img_ACR(round(centroid(2)):end,:,1),[w_point+round(33/res_ACR(2))+m w_point+round(33/res_ACR(2))+m],[n_point s_point]);
    baseline(m) = min(line_prof(:,m));
    temp_pks = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    if length(temp_pks) < 4
        pks(m,:) = zeros(1,4);
        locs(m,:) = zeros(1,4);
        w(m,:) = zeros(1,4);
        prom(m,:) = zeros(1,4);
    else
        [pks(m,:),locs(m,:),w(m,:),prom(m,:)] = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    end
end

pk_spacing = sort(locs,2);
diff_pk_spacing = diff(pk_spacing,1,2);
spacing_test = diff_pk_spacing <= ceil(hole_size_lr(2)); % peaks should be relatively close
spacing_pass = sum(spacing_test,2) >= 1; % pass if at least two peaks have an adequately spaced neighbour
w_test = 0.25*hole_size_lr(2) < w & w < 1.75*hole_size_lr(2); % peaks should have comparable width to feature
pk_test = pks > 0.4*max(line_prof(:)); % pass if peaks satisfy threshold
prom_test = prom > 0.2*max(line_prof(:)); % pass if prominences satisfy threshold

lr_resolv_10 = find(sum(spacing_pass.*w_test.*pk_test.*prom_test,2) >= 3); % output line profiles that pass
temp = lr_resolv_10; % set temp array

if ~isempty(find(diff(temp)>1))
    grouping = find(diff(temp)>1); % find clusters of line profiles that pass
    group = temp(1:grouping(1)); % find first cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    lr_resolv_10(1) = group(ind(1)); % take greatest peak from the cluster
    if length(grouping) >= 2
        for m = 1:length(grouping)-1
            group = temp(grouping(m)+1:grouping(m+1)); % find second and third clusters
            [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
            lr_resolv_10(m+1) = group(ind(1)); % take greatest peak from each cluster
        end
    end
    group = temp(grouping(end)+1:end); % find last cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    lr_resolv_10(length(grouping)+1) = group(ind(1)); % take greatest peak from the cluster

    lr_resolv_10(length(grouping)+2:end) = [];
    lr_resolv_10 = unique(lr_resolv_10);

    if length(lr_resolv_10) < 4 % check if large cluster has been overlooked
        unused = temp(~ismember(temp,lr_resolv_10)); % find passing line profiles not used
        cluster_check = sum(abs(unused - lr_resolv_10') > 1,2) == length(lr_resolv_10); % check which do not neighbour a used line profile
        [~,ind] = sort(sum(pks(unused(cluster_check),:),2),'descend'); % sort by peak height
        if ~isempty(ind)
            lr_resolv_10(end+1) = unused(ind(1)); % take greatest peak
            lr_resolv_10 = sort(lr_resolv_10);
        end
    end
elseif ~isempty(lr_resolv_10) && length(lr_resolv_10) > 4
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    lr_resolv_10 = temp(ind(1:4)); % take greatest peak
    lr_resolv_10 = sort(lr_resolv_10);
elseif ~isempty(temp) 
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    lr_resolv_10 = temp(ind(1)); % take greatest peak
end

if isrow(lr_resolv_10)
    lr_resolv_10 = lr_resolv_10';
end
%% 0.9 Horizontal UL Array

line_prof = [];
for m = 1:(s_point - n_point)/2
    line_prof(:,m) = improfile(img_ACR(round(centroid(2)):end,:,1),[w_point+round(45/res_ACR(1)) w_point + round(65/res_ACR(1))],[n_point+m n_point+m]);
    baseline(m) = min(line_prof(:,m));
    temp_pks = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    if length(temp_pks) < 4
        pks(m,:) = zeros(1,4);
        locs(m,:) = zeros(1,4);
        w(m,:) = zeros(1,4);
        prom(m,:) = zeros(1,4);
    else
        [pks(m,:),locs(m,:),w(m,:),prom(m,:)] = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    end
end

pk_spacing = sort(locs,2);
diff_pk_spacing = diff(pk_spacing,1,2);
spacing_test = diff_pk_spacing <= ceil(hole_size_ul(3)); % peaks should be relatively close
spacing_pass = sum(spacing_test,2) >= 1; % pass if at least two peaks have an adequately spaced neighbour
w_test = 0.25*hole_size_ul(3) < w & w < 1.75*hole_size_ul(3); % peaks should have comparable width to feature
pk_test = pks > 0.2*max(line_prof(:)); % pass if peaks satisfy threshold
prom_test = prom > 0.3*max(line_prof(:)); % pass if prominences satisfy threshold

ul_resolv_09 = find(sum(spacing_pass.*w_test.*pk_test.*prom_test.*w_test,2) == 4); % output line profiles that pass
temp = ul_resolv_09; % set temp array

if ~isempty(find(diff(temp)>1))
    grouping = find(diff(temp)>1); % find clusters of line profiles that pass

    group = temp(1:grouping(1)); % find first cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    ul_resolv_09(1) = group(ind(1)); % take greatest peak from the cluster
    if length(grouping) >= 2
        for m = 1:length(grouping)-1
            group = temp(grouping(m)+1:grouping(m+1)); % find second and third clusters
            [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
            ul_resolv_09(m+1) = group(ind(1)); % take greatest peak from each cluster
        end
    end
    group = temp(grouping(end)+1:end); % find last cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    ul_resolv_09(length(grouping)+1) = group(ind(1)); % take greatest peak from the cluster

    ul_resolv_09(length(grouping)+2:end) = [];
    ul_resolv_09 = unique(ul_resolv_09);

    if length(ul_resolv_09) < 4 % check if large cluster has been overlooked
        unused = temp(~ismember(temp,ul_resolv_09)); % find passing line profiles not used
        cluster_check = sum(abs(unused - ul_resolv_09') > 1,2) == length(ul_resolv_09); % check which do not neighbour a used line profile
        [~,ind] = sort(sum(pks(unused(cluster_check),:),2),'descend'); % sort by peak height
        if ~isempty(ind)
            ul_resolv_09(end+1) = unused(ind(1)); % take greatest peak
            ul_resolv_09 = sort(ul_resolv_09);
        end
    end
elseif ~isempty(ul_resolv_09) && length(ul_resolv_09) > 4
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    ul_resolv_09 = temp(ind(1:4)); % take greatest peak
    ul_resolv_09 = sort(ul_resolv_09);
elseif ~isempty(temp) 
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    ul_resolv_09 = temp(ind(1)); % take greatest peak
end

if isrow(ul_resolv_09)
    ul_resolv_09 = ul_resolv_09';
end
%% 0.9 Vertical LR Array

line_prof = [];
for m = 1:(s_point - n_point)/2
    line_prof(:,m) = improfile(img_ACR(round(centroid(2)):end,:,1),[w_point+round(60/res_ACR(2))+m w_point+round(60/res_ACR(2))+m],[n_point s_point]);
    baseline(m) = min(line_prof(:,m));
    temp_pks = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    if length(temp_pks) < 4
        pks(m,:) = zeros(1,4);
        locs(m,:) = zeros(1,4);
        w(m,:) = zeros(1,4);
        prom(m,:) = zeros(1,4);
    else
        [pks(m,:),locs(m,:),w(m,:),prom(m,:)] = findpeaks(line_prof(:,m)-baseline(m),'SortStr','descend','NPeaks',4);
    end
end

pk_spacing = sort(locs,2);
diff_pk_spacing = diff(pk_spacing,1,2);
spacing_test = diff_pk_spacing <= ceil(hole_size_lr(3)); % peaks should be relatively close
spacing_pass = sum(spacing_test,2) >= 1; % pass if at least two peaks have an adequately spaced neighbour
w_test = 0.25*hole_size_lr(3) < w & w < 1.5*hole_size_lr(3); % peaks should have comparable width to feature
pk_test = pks > 0.2*max(line_prof(:)); % pass if peaks satisfy threshold
prom_test = prom > 0.3*max(line_prof(:)); % pass if prominences satisfy threshold

lr_resolv_09 = find(sum(spacing_pass.*w_test.*pk_test.*prom_test,2) == 4); % output line profiles that pass
temp = lr_resolv_09; % set temp array

if ~isempty(find(diff(temp)>1))
    grouping = find(diff(temp)>1); % find clusters of line profiles that pass

    group = temp(1:grouping(1)); % find first cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    lr_resolv_09(1) = group(ind(1)); % take greatest peak from the cluster
    if length(grouping) >= 2
        for m = 1:length(grouping)-1
            group = temp(grouping(m)+1:grouping(m+1)); % find second and third clusters
            [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
            lr_resolv_09(m+1) = group(ind(1)); % take greatest peak from each cluster
        end
    end
    group = temp(grouping(end)+1:end); % find last cluster
    [~,ind] = sort(sum(pks(group,:),2),'descend'); % sort by peak height
    lr_resolv_09(length(grouping)+1) = group(ind(1)); % take greatest peak from the cluster

    lr_resolv_09(length(grouping)+2:end) = [];
    lr_resolv_09 = unique(lr_resolv_09);

    if length(lr_resolv_09) < 4 % check if large cluster has been overlooked
        unused = temp(~ismember(temp,lr_resolv_09)); % find passing line profiles not used
        cluster_check = sum(abs(unused - lr_resolv_09') > 1,2) == length(lr_resolv_09); % check which do not neighbour a used line profile
        [~,ind] = sort(sum(pks(unused(cluster_check),:),2),'descend'); % sort by peak height
        if ~isempty(ind)
            lr_resolv_09(end+1) = unused(ind(1)); % take greatest peak
            lr_resolv_09 = sort(lr_resolv_09);
        end
    end
elseif ~isempty(lr_resolv_09) && length(lr_resolv_09) > 4
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    lr_resolv_09 = temp(ind(1:4)); % take greatest peak
    lr_resolv_09 = sort(lr_resolv_09);
elseif ~isempty(temp) 
    [~,ind] = sort(sum(pks(temp,:),2),'descend'); % sort by peak height
    lr_resolv_09 = temp(ind(1)); % take greatest peak
end

if isrow(lr_resolv_09)
    lr_resolv_09 = lr_resolv_09';
end
%% Output
% Resolvable Line Profiles

imshow(img_insert(round(centroid(2)):end,:,1),[])
hold on
if ~isempty(ul_resolv_11)
    if length(ul_resolv_11) == 2
        plot([w_point w_point + round(15/res_ACR(1))]',[n_point+ul_resolv_11 n_point+ul_resolv_11]') % UL 1.1 line profiles
    else
        plot([w_point w_point + round(15/res_ACR(1))],[n_point+ul_resolv_11 n_point+ul_resolv_11])
    end
end
if ~isempty(lr_resolv_11)
    if length(lr_resolv_11) == 2
        plot([w_point+round(10/res_ACR(1))+lr_resolv_11 w_point+round(10/res_ACR(1))+lr_resolv_11]',[n_point s_point]') % LR 1.1 line profiles
    else
        plot([w_point+round(10/res_ACR(1))+lr_resolv_11 w_point+round(10/res_ACR(1))+lr_resolv_11],[n_point s_point]) 
    end
end
if ~isempty(ul_resolv_10)
    if length(ul_resolv_10) == 2
        plot([w_point+round(25/res_ACR(1)) w_point + round(45/res_ACR(1))]',[n_point+ul_resolv_10 n_point+ul_resolv_10]') % UL 1.0 line profiles
    else
        plot([w_point+round(25/res_ACR(1)) w_point + round(45/res_ACR(1))],[n_point+ul_resolv_10 n_point+ul_resolv_10])
    end
end
if ~isempty(lr_resolv_10)
    if length(lr_resolv_10) == 2
        plot([w_point+round(33/res_ACR(1))+lr_resolv_10 w_point+round(33/res_ACR(1))+lr_resolv_10]',[n_point s_point]') % LR 1.0 line profiles
    else
        plot([w_point+round(33/res_ACR(1))+lr_resolv_10 w_point+round(33/res_ACR(1))+lr_resolv_10],[n_point s_point]) 
    end
end
if ~isempty(ul_resolv_09)
    if length(ul_resolv_09) == 2
    plot([w_point+round(45/res_ACR(1)) w_point+round(65/res_ACR(1))]',[n_point+ul_resolv_09 n_point+ul_resolv_09]') % UL 0.9 line profiles
    else
        plot([w_point+round(45/res_ACR(1)) w_point+round(65/res_ACR(1))],[n_point+ul_resolv_09 n_point+ul_resolv_09]) 
    end
end
if ~isempty(lr_resolv_09)
    if length(lr_resolv_09) == 2
        plot([w_point+round(60/res_ACR(1))+lr_resolv_09 w_point+round(60/res_ACR(1))+lr_resolv_09]',[n_point s_point]') % LR 0.9 line profiles
    else
        plot([w_point+round(60/res_ACR(1))+lr_resolv_09 w_point+round(60/res_ACR(1))+lr_resolv_09],[n_point s_point]) 
    end
end

if isempty(ul_resolv_11) && isempty(lr_resolv_11)
    resolvable = 'No array pair was able to be fully resolved.';
elseif ~isempty(ul_resolv_09) && ~isempty(lr_resolv_09)
    resolvable = 'UL = 0.9mm, LR = 0.9mm';
elseif ~isempty(ul_resolv_09) && ~isempty(lr_resolv_10)
    resolvable = 'UL = 0.9mm, LR = 1.0mm';
elseif ~isempty(ul_resolv_10) && ~isempty(lr_resolv_09)
    resolvable = 'UL = 1.0mm, LR = 0.9mm';
elseif ~isempty(ul_resolv_09) && ~isempty(lr_resolv_11)
    resolvable = 'UL = 0.9mm, LR = 1.1mm';
elseif ~isempty(ul_resolv_11) && ~isempty(lr_resolv_09)
    resolvable = 'UL = 1.1mm, LR = 0.9mm';
elseif ~isempty(ul_resolv_10) && ~isempty(lr_resolv_10)
    resolvable = 'UL = 1.0mm, LR = 1.0mm';
elseif ~isempty(ul_resolv_10) && ~isempty(lr_resolv_11)
    resolvable = 'UL = 1.0mm, LR = 1.1mm';
elseif ~isempty(ul_resolv_11) && ~isempty(lr_resolv_10)
    resolvable = 'UL = 1.1mm, LR = 1.0mm';
elseif ~isempty(ul_resolv_11) && ~isempty(lr_resolv_11)
    resolvable = 'UL = 1.1mm, LR = 1.1mm';
elseif ~isempty(ul_resolv_09)
    resolvable = 'UL = 0.9mm, LR = N/A';
elseif ~isempty(ul_resolv_10)
    resolvable = 'UL = 1.0mm, LR = N/A';
elseif ~isempty(ul_resolv_11)
    resolvable = 'UL = 1.1mm, LR = N/A';
elseif ~isempty(lr_resolv_09)
    resolvable = 'UL = N/A, LR = 0.9mm';
elseif ~isempty(lr_resolv_10)
    resolvable = 'UL = N/A, LR = 1.0mm';
elseif ~isempty(lr_resolv_11)
    resolvable = 'UL = N/A, LR = 1.1mm';
end

title('Resolvable Line Profiles')