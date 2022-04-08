%% ACR Edge Identify
% by Yassine Azma (Apr 2022)
%
% This script takes the desired slice used for MTF measurement and identifies
% an appropriate edge feature. For slice 5 in the newer ACR phantom model, 
% edge locations are produced for both horizontal and vertical MTF 
% measurement. For an edge produced from slice 1, only one location is 
% produced.

function [x,y] = ACR_EdgeIdentify(img_insert,img_ACR,rot_ang,slicenum,res_ACR,obj_ACR)

if slicenum == 5
    thresh_img = ACR_Threshold(img_insert,res_ACR).*(img_insert > 0.1*max(img_insert(:)));
    label_map = bwlabel(thresh_img,4);
    if max(label_map(:)) > 1
        for k = 2:max(label_map(:))
            [row,col] = find(label_map==k);
            coords(k-1,:) = [mean(row) mean(col)];
        end
        x_pos = coords(:,2);
        y_pos = coords(:,1);

        middle = ACR_Centroid(img_ACR,obj_ACR);
        left_ind = [ACR_TargetIdentify(x_pos,y_pos,middle,res_ACR,rot_ang,'bottom left'),...
            ACR_TargetIdentify(x_pos,y_pos,middle,res_ACR,rot_ang,'middle left'),...
            ACR_TargetIdentify(x_pos,y_pos,middle,res_ACR,rot_ang,'top left')];
        left_side = round(mean(x_pos(left_ind))-8/res_ACR(2)); % col location
        left_y = mean([coords(left_ind(2),1),coords(left_ind(3),1)]);
        left_x = left_side + (middle(1)-left_y)*tand(rot_ang);

        top_ind = [ACR_TargetIdentify(x_pos,y_pos,middle,res_ACR,rot_ang,'top left'),...
            ACR_TargetIdentify(x_pos,y_pos,middle,res_ACR,rot_ang,'top middle'),...
            ACR_TargetIdentify(x_pos,y_pos,middle,res_ACR,rot_ang,'top right')];
        top_side = round(mean(y_pos(top_ind))-8/res_ACR(1)); % row location
        top_x = mean(x_pos(top_ind(2:end)));
        top_y = top_side - (middle(2)-top_x)*tand(rot_ang);

    end
    pe_dir = obj_ACR.getAttributeByName('InPlanePhaseEncodingDirection');
    if strcmp(pe_dir,'COL')
        x = [left_x, top_x]; %
        y = [left_y, top_y];
    else
        x = [top_x, left_x];
        y = [top_y, left_y];
    end
else
    centroid = ACR_Centroid(img_ACR,obj_ACR);
    c = improfile(img_insert,[centroid(2) centroid(2)],centroid(1)+[-20 20]);
    diff_c = abs(diff(c));
    [~,pk_loc] = findpeaks(diff_c,'NPeaks',1,'Threshold',0.05*max(diff_c));

    x = centroid(2);
    y = centroid(1)+(pk_loc+1)-20;
end


