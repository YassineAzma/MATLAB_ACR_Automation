%% ACR Target Identify
% by Yassine Azma (Apr 2022)
%
% This script takes slice 5 from the newer ACR model and identifies the
% point targets within it needed for finding the edge.

function ind = ACR_TargetIdentify(x_pos,y_pos,centroid,res_ACR,rot_ang,target)

dist_x = x_pos - centroid(1);
dist_y = y_pos - centroid(2);

half_t2t_dist = round(30*cosd(rot_ang)/res_ACR(1)); % half the target to target distance
left_targets = dist_x < -half_t2t_dist;
right_targets = dist_x > half_t2t_dist;
bottom_targets = dist_y > half_t2t_dist;
top_targets = dist_y < -half_t2t_dist;
middle_column_targets = ~left_targets.*~right_targets;
middle_row_targets = ~bottom_targets.*~top_targets;

switch target
    case 'bottom left'
        ind = find(bottom_targets.*left_targets);
    case 'middle left'
        ind = find(middle_row_targets.*left_targets);
    case 'top left'
        ind = find(top_targets.*left_targets);
    case 'bottom middle'
        ind = find(middle_column_targets.*bottom_targets);
    case 'middle middle'
        ind = find(middle_column_targets.*middle_row_targets);
    case 'top middle'
        ind = find(top_targets.*middle_column_targets);
    case 'bottom right'
        ind = find(bottom_targets.*right_targets);
    case 'middle right'
        ind = find(middle_row_targets.*right_targets);
    case 'top right'
        ind = find(top_targets.*right_targets);
end
