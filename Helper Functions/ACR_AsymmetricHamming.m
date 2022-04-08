%% ACR Asymmetric Hamming
% by Yassine Azma (Apr 2022)
%
% This script creates an asymmetric Hamming window with a non-zero plateau.

function hamming_window = ACR_AsymmetricHamming(N,Centre)

hamming_window = 0.54 + 0.46*cos(([1:N]-Centre)*2*pi/N);
min_val = min(hamming_window);

turning_pt = find(min_val==hamming_window);
if turning_pt < Centre
    hamming_window(1:turning_pt-1) = min_val;
else
    hamming_window(turning_pt+1:end) = min_val;
end