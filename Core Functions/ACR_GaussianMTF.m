% function eff_res = ACR_MTF(img_ACR,obj_ACR)
% close all

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_insert = squeeze(double(img_ACR(:,:,1,1))); % if yes, only process the first
    waitfor(msgbox('4D array detected. Only processing first axial series.'));
else
    img_insert = double(img_ACR(:,:,1));
end

rot_ang = ACR_FindRotation(img_ACR,obj_ACR);

if rot_ang == 0
    waitfor(msgbox('The phantom has no detectable tilt so MTF results will not be accurate. Click OK to continue.'));
end

level = 0.5;
res_ACR = ACR_RetrievePixelSpacing(obj_ACR);

imshow(img_insert,[],'InitialMagnification',300) % magnify
roi = createMask(drawrectangle); % create rectangular ROI

%% 2D Gaussian Fitting
[row,col] = find(roi);

min_x = min(col);
min_y = min(row);
max_x = max(col);
max_y = max(row);

roi_crop = imcrop(img_insert,[min_x min_y max_x-min_x max_y-min_y]); % crop image
if size(roi_crop,2) > size(roi_crop,1) % Find edge axis based on which ROI dimension is longer
    roi_crop = roi_crop';
    trans = 1; % 
end

[X,Y] = meshgrid(1:size(roi_crop,2),1:size(roi_crop,1));

thresh_roi_crop = roi_crop > 0.8*max(roi_crop(:));
vh = mean(roi_crop(thresh_roi_crop));
vl = min(roi_crop(~thresh_roi_crop));

fun = @(a,X,Y) (vh-vl)*normcdf(Y,a(1)+a(2)*X,0.5)+vl;
obj_fun = @(a) sum((roi_crop - fun(a,X,Y)).^2,'all');
[par,fit] = fminsearch(@(a) obj_fun(a),[0 0]);

%% Resampling Along Edge
edge_slope = -1/par(2);

if trans == 1
    resample_roi_crop = imresize(roi_crop,[size(roi_crop,1),8*size(roi_crop,2)]);
else
    resample_roi_crop = imresize(roi_crop,[8*size(roi_crop,1),size(roi_crop,2)]);
end

mid_loc = [ceil(size(resample_roi_crop,1)/2) ceil(size(resample_roi_crop,2)/2)];
[X_resample,Y_resample] = meshgrid(1:size(resample_roi_crop,2),1:size(resample_roi_crop,1));

diffY = flipud(Y_resample-1);
X_prime = X_resample + diffY*edge_slope;

x_range = [floor(min(X_prime,[],'all')), ceil(max(X_prime,[],'all'))];

erf = [];
for k = x_range(1):x_range(2)
    index = k+abs(x_range(1))+1;
    if k < x_range(2)
        erf(index) = mean(resample_roi_crop((X_prime > k) .* (X_prime <= k+1)==1));
    end
end

%% LSF
hamming_filt = hamming(length(erf)-1)'; % hamming filter
lsf = diff(erf); % calculate LSF
lsf = hamming_filt.*lsf;

%% MTF
N = size(lsf,2); % length of LSF

if mod(N,2)==0
    n = -N/2:N/2-1; % even
else
    n = -(N-1)/2:(N-1)/2; % odd
end

Fs = 1/(rms(res_ACR)*(1/8)); % sampling spatial frequency (lp/mm)
T = N/Fs; % 
freq = n/T; % convert spatial period to frequency  
MTF = abs(fftshift(fft(lsf))); % perform FFT for MTF
MTF = MTF./max(MTF);
zero_freq = find(freq==0); % find zero frequency

figure
subplot(2,2,1)
plot(erf,'b.')
% Label axes
xlabel( 'Pixel', 'Interpreter', 'none' );
ylabel( 'Intensity', 'Interpreter', 'none' );
grid on
title('Edge Response Function');

subplot(2,2,2)
plot(lsf,'b.')
% Label axes
xlabel( 'Pixel', 'Interpreter', 'none' );
ylabel( 'Intensity', 'Interpreter', 'none' );
grid on
title('Line Spread Function');

subplot(2,2,3)
plot(freq,MTF,'b.')
% Label axes
ylim([0 1])
xlabel( 'Spatial Frequency (lp/mm)', 'Interpreter', 'none' );
ylabel( '', 'Interpreter', 'none' );
grid on
title('Modulation Transfer Function');
legend('Gaussian Fitting Data', 'Location', 'NorthEast', 'Interpreter', 'none' );

MTF_corr = MTF(zero_freq:end); % discard symmetrical negative freqs of MTF
freq_corr = freq(zero_freq:end); % same as above for spatial frequencies

freq_interp = interp1(MTF_corr,freq_corr,0:0.005:1,'pchip'); % interpolate MTF
eq_lp = freq_interp(0:0.005:1==level); % Closest spatial frequency to a MTF of 0.5

eff_res = (1/(eq_lp*2)); % output effective resolution for the line profile

subplot(2,2,4)
plot(freq_corr,MTF_corr,'b.')
% Label axes
xlim([0 5])
ylim([0 1])
xlabel( 'Spatial Frequency (lp/mm)', 'Interpreter', 'none' );
ylabel( '', 'Interpreter', 'none' );
grid on
title('Modulation Transfer Function');
legend(['Original Data - ' num2str(round(eff_res,2)) 'mm @ ' num2str(100*level) '%'], ...
    'Location', 'NorthEast', 'Interpreter', 'none' );