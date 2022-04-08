% by Yassine Azma (Mar 2022)
%
% This script takes the ACR dataset, identifies the ACR model and uses the
% most appropriate slice for an MTF measurement. Alternatively, a
% rectangular ROI can be chosen by the user. Using a 2D surface fit, 
% the edge slope is determined in a noise-robust manner. Using this slope, 
% line profiles are acquired parallel to the edge to obtain the edge
% response function. The line spread function is determined and smoothed
% with a Hamming filter before a Fourier Transform to obtain the MTF.

function eff_res = ACR_GaussianMTF(img_ACR,obj_ACR,mode,options)

model = ACR_IdentifyModel(img_ACR);

if strcmp(model,'old')
    slicenum = 1;
else
    slicenum = 5;
end

if size(img_ACR,4) > 1 % check if input array contains multiple ACR series
    img_insert = squeeze(double(img_ACR(:,:,slicenum,1))); % if yes, only process the first
    waitfor(msgbox('4D array detected. Only processing first axial series.'));
else
    img_insert = double(img_ACR(:,:,slicenum));
end

rot_ang = ACR_FindRotation(img_ACR,obj_ACR);

if abs(rot_ang) < 1
    dlgopts.Interpreter = 'tex';
    dlgopts.Default = 'Yes';
    answer = questdlg(['The phantom has no significant tilt (|\theta| < 1' char(176) ') ' ...
        '. MTF results will not be accurate. Do you wish to continue?'],'MTF Accuracy Warning',...
        'Yes','No',dlgopts);
    switch answer
        case 'Yes'
        case 'No'
            eff_res = [];
            return
    end
end

level = 0.5; % MTF ratio to report effective resolution - MTF50 is a common standard
res_ACR = ACR_RetrievePixelSpacing(obj_ACR);

if strcmp(mode,'auto')
    [x,y] = ACR_EdgeIdentify(img_insert,img_ACR,rot_ang,slicenum,res_ACR,obj_ACR);
else
    figure
    imshow(img_insert,[],'InitialMagnification',300) % magnify
    colormap('parula')
    [x,y] = ginput(1);
    close
end

for q = 1:length(x)
    if rot_ang < 1
        width = round(17*size(img_insert,1)/256)
    else
        width = round(17*size(img_insert,1)/256*(5/rot_ang));
    end
    coords = [x(q)-width/2 y(q)-width/2];
    roi_crop = imcrop(img_insert,[coords(1) coords(2) width width]); % crop image
    if strcmp(options.ModelSensitivity,'yes')
        sens_crop = imcrop(ACR_SensitivityMap(img_ACR),[coords(1) coords(2) width width]); % crop sensitivity map
    end
    
    % 2D Gaussian Fitting
    edge_sum_cols = sum(roi_crop,1);
    edge_sum_rows = sum(roi_crop,2);
    pk_cols = findpeaks(abs(diff(edge_sum_cols)),'NPeaks',1,'SortStr','descend');
    pk_rows = findpeaks(abs(diff(edge_sum_rows)),'NPeaks',1,'SortStr','descend');

    thresh_roi_crop = roi_crop > 0.6*max(roi_crop(:));

    if pk_rows > pk_cols % Edge is vertical
        edge_type = 'vertical';
        edge_dir = sum(thresh_roi_crop,1); % direction edge increases, increasing means downward

        naive_lsf = abs(diff(sum(thresh_roi_crop,2)))>1;
        edge_test = diff(find(naive_lsf==0));
        edge_begin = find(edge_test > 1)+1;
        edge_loc = [edge_begin-1 edge_begin + edge_test(edge_begin-1)-1];
        if edge_dir(end) > edge_dir(1)
            dir = 'downward';
        else
            dir = 'upward';
        end
    else
        edge_type = 'horizontal';
        edge_dir = sum(thresh_roi_crop,2); % direction edge increases, increasing means leftward

        naive_lsf = abs(diff(sum(thresh_roi_crop,1)))>1;
        edge_test = diff(find(naive_lsf==0));
        edge_begin = find(edge_test > 1)+1;
        edge_loc = [edge_begin-1 edge_begin + edge_test(edge_begin-1)-1];
        if edge_dir(end) > edge_dir(1)
            dir = 'leftward';
        else
            dir = 'rightward';
        end
    end

    [X,Y] = meshgrid(1:size(roi_crop,2),1:size(roi_crop,1));

    vh = max(roi_crop(thresh_roi_crop));
    vl = 20+min(roi_crop(~thresh_roi_crop));

    switch [dir ' ' edge_type]
        case 'downward vertical'
            ft = fittype( '(h-l)*normcdf(x,a+b*y,0.5)+l', 'independent', {'x', 'y'}, 'dependent', 'z' );
        case 'upward vertical'
            ft = fittype( '(h-l)*normcdf(-x,a+b*y,0.5)+l', 'independent', {'x', 'y'}, 'dependent', 'z' );
        case 'leftward horizontal'
            ft = fittype( '(h-l)*normcdf(y,a+b*x,0.5)+l', 'independent', {'x', 'y'}, 'dependent', 'z' );
        case 'rightward horizontal'
            ft = fittype( '(h-l)*normcdf(-y,a+b*x,0.5)+l', 'independent', {'x', 'y'}, 'dependent', 'z' );
    end

    [xData, yData, zData] = prepareSurfaceData( X, Y, roi_crop );

    % Set up fittype and options.
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf vh vl];
    opts.StartPoint = [0 0 vh vl];
    opts.Upper = [Inf Inf vh vl];

    % Fit model to data.
    [fitresult, gof] = fit( [xData, yData], zData, ft, opts );

    edge_dir = sum(thresh_roi_crop,2);
    par = coeffvalues(fitresult);

    if edge_dir(end) > edge_dir(1)
        edge_slope = 1/par(2);
    else
        edge_slope = -1/par(2);
    end

    % Resampling Along Edge
    resamp_factor = 8;
    if strcmp(edge_type,'horizontal')
        resample_roi_crop = imresize(roi_crop,[size(roi_crop,1),resamp_factor*size(roi_crop,2)]);
        if strcmp(options.ModelSensitivity,'yes')
            resample_sens_crop = imresize(sens_crop,[size(roi_crop,1),resamp_factor*size(roi_crop,2)]);
        end
    else
        resample_roi_crop = imresize(roi_crop,[resamp_factor*size(roi_crop,1),size(roi_crop,2)]);
        if strcmp(options.ModelSensitivity,'yes')
            resample_sens_crop = imresize(sens_crop,[resamp_factor*size(roi_crop,1),size(roi_crop,2)]);
        end
    end

    mid_loc = [ceil(size(resample_roi_crop,1)/2) ceil(size(resample_roi_crop,2)/2)];
    [X_resample,Y_resample] = meshgrid(1:size(resample_roi_crop,2),1:size(resample_roi_crop,1));

    if strcmp(edge_type,'horizontal')
        diffY = flipud(Y_resample-1)-mid_loc(1);
        X_prime = X_resample + resamp_factor*diffY*edge_slope;

        x_range = [floor(min(X_prime,[],'all')), ceil(max(X_prime,[],'all'))];

%             figure
        erf = [];
        for k = x_range(1):x_range(2)
            index = k+abs(x_range(1))+1;
            if k < x_range(2)
                erf(index) = mean(resample_roi_crop((X_prime >= k) .* (X_prime < k+1)==1));
                if strcmp(options.ModelSensitivity,'yes')
                    sens(index) = mean(resample_sens_crop((X_prime >= k) .* (X_prime < k+1)==1));
                end
                n_inside_roi(index) = nnz(resample_roi_crop((X_prime >= k) .* (X_prime < k+1)==1));
            end
%                     imshowpair(X_prime>=k&X_prime<k+1,resample_roi_crop)
%                     drawnow
        end
%         close
    else
        diffX = (size(X_resample,1)-X_resample)-mid_loc(2);
        Y_prime = flipud(Y_resample)+edge_slope*resamp_factor.*diffX;

        y_range = [floor(min(Y_prime,[],'all')), ceil(max(Y_prime,[],'all'))];

%             figure
        erf = [];
        for k = y_range(1):y_range(2)
            index = k+abs(y_range(1))+1;
            if k < y_range(2)
                erf(index) = mean(resample_roi_crop((Y_prime >= k) .* (Y_prime < k+1)==1));
                if strcmp(options.ModelSensitivity,'yes')
                    sens(index) = mean(resample_sens_crop((Y_prime >= k) .* (Y_prime < k+1)==1));
                end
                n_inside_roi(index) = nnz(resample_roi_crop((Y_prime >= k) .* (Y_prime < k+1)==1)); % number of non-zero values in line profile
            end
%                     imshowpair(Y_prime>=k&Y_prime<k+1,resample_roi_crop)
%                     drawnow
        end
%         close
    end
    % FIT SIGMOID TO ERF
    erf = erf(n_inside_roi == max(n_inside_roi)); % remove truncated lines (i.e. lines that would continue outside cropped roi)

    [~,pk_centre] = findpeaks(abs(diff(erf)),'SortStr','descend','NPeaks',1);
    [~,locs] = findpeaks(abs(diff(erf)));
%     inflect_pts = find(diff(sign(diff(erf)))) + 1;
%     left_side = inflect_pts(find(inflect_pts < pk_centre,1,'last'));
    [~,closestIndex] = min(abs(locs(locs < pk_centre)-pk_centre));
    left_side = locs(closestIndex);
    if isempty(left_side)
        left_side = 1;
    end
%     right_side = inflect_pts(find(inflect_pts > pk_centre,1,'first'));
    sub_pks_right = locs > pk_centre;
    [~,closestIndex] = min(abs(locs(sub_pks_right)-pk_centre));
    right_side = locs(find(sub_pks_right,1)+closestIndex-1);
    if isempty(right_side)
        right_side = length(erf);
    end
    if strcmp(options.ModelSensitivity,'yes')
        sens = sens(n_inside_roi == max(n_inside_roi));
        sens = sens./mean(sens(left_side:right_side));
        weights = ones(1,length(erf));
    else
        weights = [0.5*ones(1,left_side) ones(1,right_side-left_side) 0.5*ones(1,length(erf)-right_side)];
    end
    pts = 1:length(erf);

    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';

    % Set up fittype and options.
    if strcmp(options.ModelSensitivity,'yes')
        [xData, yData, zData, weightsData] = prepareSurfaceData(pts, sens, erf, weights);

        ft = fittype( 'a + b/(1+exp(c*(x-d)))^e + f*y*x', 'independent', {'x', 'y'}, 'dependent', 'z' );
        opts.Lower = [-Inf -Inf -Inf 0 1 -Inf];
        opts.StartPoint = [0 max(erf) 0.2 pk_centre 1 0];
        opts.Weights = weightsData;

        [fitresult_erf, ~] = fit( [xData, yData], zData, ft, opts );
        erf_sig = feval(fitresult_erf,xData,yData);
    else
        [xData, yData, weightsData] = prepareCurveData(pts, erf, weights);

        ft = fittype( 'a + b/(1+exp(c*(x-d)))^e', 'independent', 'x', 'dependent', 'y' );
        opts.Lower = [-Inf -Inf -Inf 0 1];
        opts.StartPoint = [0 max(erf) 0.2 pk_centre 1];
        opts.Weights = weightsData;

        % Fit model to data.
        [fitresult_erf, ~] = fit( xData, yData, ft, opts );
        erf_sig = feval(fitresult_erf,xData);
    end
    
    %% LSF
    lsf = diff(erf); % calculate LSF
    N = size(lsf,2); % length of LSF
    hamming_filt = ACR_AsymmetricHamming(N,pk_centre+1);
    lsf = hamming_filt.*lsf;

    lsf_sig = diff(erf_sig); % calculate LSF
    lsf_sig = hamming_filt.*lsf_sig';

    weights_lsf = weights(2:end);

    % FIT GAUSSIAN TO LSF
    [xData, yData, weightsData] = prepareCurveData( [], lsf, weights_lsf );
    ft = fittype( 'a*exp(-((x-b)/c)^2)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf 0 -Inf];
    opts.Robust = 'Bisquare';
    if abs(min(lsf)) < max(lsf)
        opts.StartPoint = [max(lsf) pk_centre 4];
    else
        opts.StartPoint = [min(lsf) pk_centre 4];
    end
    opts.Weights = weightsData;

    [fitresult_lsf, ~] = fit( xData, yData, ft, opts );
    lsf_gauss = feval(fitresult_lsf,xData);

    % MTF
    if mod(N,2)==0
        n = -N/2:N/2-1; % even
    else
        n = -(N-1)/2:(N-1)/2; % odd
    end

    Fs = 1/(rms(res_ACR)*(1/resamp_factor)); % sampling spatial frequency (lp/mm)
    freq = n*Fs/N; % convert spatial period to frequency
    MTF = abs(fftshift(fft(lsf))); % perform FFT for MTF
    MTF = MTF./max(MTF);
    MTF_sig = abs(fftshift(fft(lsf_sig))); % perform FFT for MTF
    MTF_sig = MTF_sig./max(MTF_sig);
    MTF_gauss = abs(fftshift(fft(lsf_gauss)));
    MTF_gauss = MTF_gauss./max(MTF_gauss);

    zero_freq = find(freq==0); % find zero frequency

    MTF_corr = MTF(zero_freq:end); % discard symmetrical negative freqs of MTF
    MTF_sig_corr = MTF_sig(zero_freq:end);
    MTF_gauss_corr = MTF_gauss(zero_freq:end);
    freq_corr = freq(zero_freq:end); % same as above for spatial frequencies

    freq_interp = interp1(MTF_corr,freq_corr,0:0.005:1,'pchip'); % interpolate MTF
    freq_interp_sig = interp1(MTF_sig_corr,freq_corr,0:0.005:1,'pchip'); % interpolate MTF
    [~,ind] = unique(MTF_gauss_corr);
    freq_interp_gauss = interp1(MTF_gauss_corr(ind),freq_corr(ind),0:0.005:1,'pchip');
    eq_lp = freq_interp(0:0.005:1==level); % Closest spatial frequency to a MTF of 0.5
    eq_lp_sig = freq_interp_sig(0:0.005:1==level);
    eq_lp_gauss = freq_interp_gauss(0:0.005:1==level);
    eff_res(q) = (1/(eq_lp*2)); % output effective resolution for the line profile
    eff_res_sig = (1/(eq_lp_sig*2));
    eff_res_gauss = (1/(eq_lp_gauss*2));

    if strcmp(options.SuppressFigures,'no')
        figure
        subplot(2,2,1)
        imagesc(img_insert)
        axis off
        axis image
        colormap('gray')
        hold on
        rectangle('Position',[coords(1) coords(2) width width],'EdgeColor','r')
        hold on
        switch [dir ' ' edge_type]
            case 'downward vertical'
                plot(coords(1)+[0:width],coords(2)+mean(edge_loc)-edge_slope*[0:width],'y--','LineWidth',1.2)
            case 'upward vertical'
                plot(coords(1)+[0:width],coords(2)+mean(edge_loc)-edge_slope*[0:width],'y--','LineWidth',1.2)
            case 'leftward horizontal'
                plot(coords(1)+mean(edge_loc)+edge_slope*[0:width],coords(2)+[0:width],'y--','LineWidth',1.2)
            case 'rightward horizontal'
                plot(coords(1)+mean(edge_loc)+edge_slope*[0:width],coords(2)+[0:width],'y--','LineWidth',1.2)
        end

        subplot(2,2,2)
        plot(erf,'b.','MarkerSize',8)
        hold on
        plot(erf_sig,'r-')
        % Label axes
        xlabel( 'Pixel', 'Interpreter', 'none' );
        ylabel( 'Intensity', 'Interpreter', 'none' );
        grid on
        legend('Raw','Weighted Sigmoid Fit of ERF','Location','best')
        title('Edge Response Function');

        subplot(2,2,3)
        plot(lsf,'b.','MarkerSize',8)
        hold on
        plot(lsf_sig,'r-')
        hold on
        plot(lsf_gauss,'k-')
        % Label axes
        xlabel( 'Pixel', 'Interpreter', 'none' );
        ylabel( 'Intensity', 'Interpreter', 'none' );
        grid on
        legend('Raw','Weighted Sigmoid Fit of ERF','Weighted Gaussian Fit of LSF','Location','best')
        title('Line Spread Function');

        subplot(2,2,4)
        plot(freq_corr,MTF_corr,'b.','MarkerSize',8)
        hold on
        plot(freq_corr,MTF_sig_corr,'r-')
        hold on
        plot(freq_corr,MTF_gauss_corr,'k-')
        % Label axes
        ylim([0 1])
        xlim([0 2/(2*rms(res_ACR))])
        xlabel( 'Spatial Frequency (lp/mm)', 'Interpreter', 'none' );
        ylabel( '', 'Interpreter', 'none' );
        grid on
        title('Modulation Transfer Function');
        legend(['Raw Data - ' num2str(round(eff_res(q),2)) 'mm @ ' num2str(100*level) '%'], ...
            ['Weighted Sigmoid Fit of ERF - ' num2str(round(eff_res_sig,2)) 'mm @ ' num2str(100*level) '%'], ...
            ['Weighted Gaussian Fit of LSF - ' num2str(round(eff_res_gauss,2)) 'mm @ ' num2str(100*level) '%'], ...
            'Location', 'NorthEast', 'Interpreter', 'none' );
    end
end