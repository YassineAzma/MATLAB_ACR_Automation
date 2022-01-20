%% ACR MTF
% by Yassine Azma (Dec 2021)
% 
% This script allows...
% Run ACR_DataSort before running this. Delete any non-image files from the
% folder - need to include some form of file filtering!
%% ROI Selection
% function eff_res = ACR_MTF(img_ACR,obj_ACR,level)

level = 0.5;
if isempty(obj_ACR.getAttributeByName('PixelSpacing')) % Multi-frame check
    list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');
    res_ACR = list.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
else
    res_ACR = obj_ACR.getAttributeByName('PixelSpacing'); % retrieve ACR in-plane resolution
end

img_insert = img_ACR(:,:,1); % resolution insert image

imshow(img_insert,[],'InitialMagnification',300) % magnify
roi = createMask(drawrectangle); % create rectangular ROI
hold on

%% Edge Response Functions

% 
vals = roi.*img_insert; % apply ROI mask
 
[row,col] = find(vals); % index roi

min_row = min(row); %min row (y1)
max_row = max(row); %max row (y2)

min_col = min(col); %min col (x1)
max_col = max(col); %max col (x2)

interp_factor = 10;

v_lines = vals(min_row:max_row,min_col:max_col); % vertical line array
h_lines = vals(min_row:max_row,min_col:max_col)'; % horizontal line array

lines_interp = [];

if mean(max(v_lines)-min(v_lines)) > mean(max(h_lines)-min(h_lines))
    for k = 1:size(v_lines,2)
        lines_interp(:,k) = interp1(v_lines(:,k),1:(1/interp_factor):length(v_lines(:,k)),'pchip'); %piece-wise shape preserving interpolation
    end
else
    for k = 1:size(h_lines,2)
        lines_interp(:,k) = interp1(h_lines(:,k),1:(1/interp_factor):length(h_lines),'pchip');
    end
end

%% EDGE ANGLE AND ANGULATE LINE PROFILES
changept = zeros(1,size(lines_interp,2));
for k = 1:size(lines_interp,2)
    temp = findchangepts(lines_interp(:,k),'MaxNumChanges',1);
    if isempty(temp)
        changept(k) = NaN;
    else
        changept(k) = temp;
    end
end
% changept = rmoutliers(changept);

edge_angle = atan((max(changept)-min(changept))*(1/interp_factor)/size(lines_interp,2));
if changept(1) > changept(end)
    edge_angle = -edge_angle;
end

rot_matrix = [cos(edge_angle) -sin(edge_angle); sin(edge_angle) cos(edge_angle)];

coords = [];
rot_coords = [];
rot_lines = [];
rot_lines_interp = [];

if mean(max(v_lines)-min(v_lines)) > mean(max(h_lines)-min(h_lines))
    line_centre = min(row) + (max(row)-min(row))/2;
    for k = 1:length(min_col:max_col)+1
        coords(1,:,k) = repmat(min_col+(k-1),1,length(min_row:max_row));
        coords(2,:,k) = [min_row:max_row];
        rot_coords(:,:,k) = rot_matrix*coords(:,:,k);
        rot_coords(1,:,k) = rot_coords(1,:,k);
        rot_coords(2,:,k) = rot_coords(2,:,k);

        rot_lines(:,k) = improfile(img_insert,rot_coords(1,:,k),rot_coords(2,:,k));
    end
else
    line_centre = min(col) + (max(col)-min(col))/2;
    for k = 1:length(min_row:max_row)
        coords(1,:,k) = [min_col:max_col];
        coords(2,:,k) = repmat(min_row+(k-1),1,length(min_col:max_col));
        rot_coords(:,:,k) = rot_matrix*coords(:,:,k);
        rot_coords(1,:,k) = rot_coords(1,:,k);
        rot_coords(2,:,k) = rot_coords(2,:,k);

        rot_lines(:,k) = improfile(img_insert,rot_coords(1,:,k),rot_coords(2,:,k));
    end
end

for k = 1:size(v_lines,2)
    rot_lines_interp(:,k) = interp1(rot_lines(:,k),1:(1/interp_factor):length(rot_lines(:,k)),'pchip'); %piece-wise shape preserving interpolation
end

for k = 1:size(rot_coords,3)
    plot(rot_coords(1,:,k),rot_coords(2,:,k),'b')
    hold on
end
hold off
%% LINE PROFILES PERPENDICULAR TO EDGE

% rot_img = imrotate(img_insert,-rad2deg(edge_angle),'bicubic','crop');
% vals = roi.*rot_img; % apply ROI mask
% [row,col] = find(vals); % index roi
% 
% min_row = min(row); %min row (y1)
% max_row = max(row); %max row (y2)
% 
% min_col = min(col); %min col (x1)
% max_col = max(col); %max col (x2)
% 
% v_lines = vals(min_row:max_row,min_col:max_col); % vertical line array
% h_lines = vals(min_row:max_row,min_col:max_col)'; % horizontal line array
% 
% v_lines_interp = [];
% h_lines_interp = [];

% for k = 1:size(v_lines,2)
%     lines_interp(:,k) = interp1(v_lines(:,k),1:(1/interp_factor):length(v_lines(:,k)),'pchip'); %piece-wise shape preserving interpolation
% end

lines_interp = sum(rot_lines_interp,2)/size(rot_lines_interp,2);
lines_interp = lines_interp - min(lines_interp);

inflect_pts = findchangepts(lines_interp,'MaxNumChanges',2);
weights = [0.5*ones(1,inflect_pts(1)-1) ones(1,diff(inflect_pts)+1) 0.5*ones(1,length(lines_interp)-inflect_pts(2))];
% FIT SIGMOID TO ERF
x = 1:length(lines_interp);
[xData, yData, weightsData] = prepareCurveData(x, lines_interp, weights);

% Set up fittype and options.
ft = fittype( 'a + b/(1+exp(c*(x-d)))^e', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf 0 0];
opts.Robust = 'Bisquare';
opts.StartPoint = [0 mean(lines_interp) 0.2 findchangepts(lines_interp,'MaxNumChanges',1) 1];
opts.Weights = weightsData;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
ERF = feval(fitresult,xData);

figure
% Plot fit with data.
subplot(2,2,1)
h = plot( fitresult, xData, yData );
legend( h, 'Original Data', 'Weighted Sigmoid Fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Pixel', 'Interpreter', 'none' );
ylabel( 'Intensity', 'Interpreter', 'none' );
grid on
title('Edge Response Function');

% CREATE LSF
hamming_filt = hamming(length(lines_interp)-1); % hamming filter

LSF = (diff(lines_interp))/(res_ACR(2)/interp_factor); % calculate LSF
LSF = hamming_filt.*LSF;
LSF = LSF./sum(LSF,1); % normalize

LSF_sig = (diff(ERF))/(res_ACR(2)/interp_factor); % calculate LSF
LSF_sig = hamming_filt.*LSF_sig;
LSF_sig = LSF_sig./sum(LSF_sig,1); % normalize

subplot(2,2,2)
plot(LSF,'b.')
hold on
plot(LSF_sig,'r-')
% Label axes
xlabel( 'Pixel', 'Interpreter', 'none' );
ylabel( 'Intensity', 'Interpreter', 'none' );
grid on
title('Line Spread Function');
legend('Original Data', 'Weighted Sigmoid Fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
%% Modulation Transfer Function
 
N = size(LSF,1); % length of LSF

if mod(N,2)==0
    n = -N/2:N/2-1; % even
else
    n = -(N-1)/2:(N-1)/2; % odd
end

Fs = 1/(rms(res_ACR)*(1/interp_factor)); % sampling spatial frequency (lp/mm)
T = N/Fs; % 
freq = n/T; % convert spatial period to frequency  
MTF = abs(fftshift(fft(LSF))); % perform FFT for MTF
MTF_sig = abs(fftshift(fft(LSF_sig))); % perform FFT for MTF
zero_freq = find(freq==0); % find zero frequency

subplot(2,2,3)
plot(freq,MTF,'b.')
hold on
plot(freq,MTF_sig,'r-')
% Label axes
ylim([0 1])
xlabel( 'Spatial Frequency (lp/mm)', 'Interpreter', 'none' );
ylabel( '', 'Interpreter', 'none' );
grid on
title('Modulation Transfer Function');
legend('Original Data', 'Weighted Sigmoid Fit', 'Location', 'NorthEast', 'Interpreter', 'none' );

MTF_corr = MTF(zero_freq:end); % discard symmetrical negative freqs of MTF
MTF_corr_sig = MTF_sig(zero_freq:end); % discard symmetrical negative freqs of MTF
freq_corr = freq(zero_freq:end); % same as above for spatial frequencies

freq_interp = [];
freq_interp = interp1(MTF_corr,freq_corr,0:0.005:1,'pchip'); % interpolate MTF
freq_interp_sig = interp1(MTF_corr_sig,freq_corr,0:0.005:1,'pchip'); % interpolate MTF
eq_lp = freq_interp(0:0.005:1==level); % Closest spatial frequency to a MTF of 0.5
eq_lp_sig = freq_interp_sig(0:0.005:1==level); % Closest spatial frequency to a MTF of 0.5

eff_res = (1/(eq_lp*2)); % output effective resolution for the line profile
eff_res_sig = (1/(eq_lp_sig*2)); % output effective resolution for the sigmoid ERF

subplot(2,2,4)
plot(freq_corr,MTF_corr,'b.')
hold on
plot(freq_corr,MTF_corr_sig,'r-')
% Label axes
xlim([0 5])
ylim([0 1])
xlabel( 'Spatial Frequency (lp/mm)', 'Interpreter', 'none' );
ylabel( '', 'Interpreter', 'none' );
grid on
title('Modulation Transfer Function');
legend(['Original Data - ' num2str(round(eff_res,2)) 'mm @ ' num2str(100*level) '%'], ...
    ['Weighted Sigmoid Fit - ' num2str(round(eff_res_sig,2)) 'mm @ ' num2str(100*level) '%'], ...
    'Location', 'NorthEast', 'Interpreter', 'none' );
