%% ACR FWHM
% by Yassine Azma (Jan 2022)
%
% This script uses the line profiles of the ramps within the insert to 
% determine the full-width half-maximum of the respective profiles.

function [pk,pk_index,indices] = ACR_FWHM(line_prof)

[pk, pk_index] = maxk(line_prof,1);
boundaries = findchangepts(line_prof,'MaxNumChanges',2);
diff_val = abs(line_prof - 0.5*pk);

if length(boundaries)==1
    % identify side
    if boundaries - pk_index > 0
        [~,r_ind] = sort(diff_val(pk_index:end));
        bound_rdist = abs(pk_index+r_ind(1:5)-boundaries); % closest 2 indices to half-maximum
        [~,rs_ind] = sort(bound_rdist);
    
        right_side = pk_index+r_ind(rs_ind(1));

        left_side = 2*pk_index-right_side;
    else
        [~,l_ind] = sort(diff_val(1:pk_index));
        bound_ldist = abs(l_ind(1:5)-boundaries); % closest 2 indices to half-maximum
        [~,ls_ind] = sort(bound_ldist);
    
        left_side = l_ind(ls_ind(1));

        right_side = 2*pk_index-left_side;
    end
else
    % left side
    [~,l_ind] = sort(diff_val(1:pk_index));
    del_ind = l_ind < boundaries(1);
    l_ind(del_ind) = [];
    bound_ldist = abs(l_ind(1:5)-boundaries(1)); % closest 2 indices to half-maximum
    [~,ls_ind] = sort(bound_ldist);
    
    left_side = l_ind(ls_ind(1));
    
    % right side
    [~,r_ind] = sort(diff_val(pk_index:end));
    del_ind = r_ind > boundaries(2);
    r_ind(del_ind) = [];
    bound_rdist = abs(pk_index+r_ind(1:5)-boundaries(2)); % closest 2 indices to half-maximum
    [~,rs_ind] = sort(bound_rdist);
    
    right_side = pk_index+r_ind(rs_ind(1));
end

indices = [left_side right_side];
