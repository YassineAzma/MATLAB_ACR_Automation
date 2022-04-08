%% ACR Data Sort
% by Yassine Azma (Oct 2021)
% 
% This script allows the user to select the two separate folders containing
% the sagittal localiser and the axial series respectively. The images are
% then reordered, in the event that an interleaved sequence was selected,
% and the DICOM metadata for later extraction of the pixel spacing are produced. 

function [img_loc, img_ACR, obj_loc, obj_ACR] = ACR_DataSort(options)

%% Localiser
if strcmp(options.IncludeLocaliser,'yes')
    switch options.Orientation
        case 'axial'
            sag_dir_cosines = [0 1 0 0 0 -1]'; % direction cosines for axial positioning
        case 'coronal'
            sag_dir_cosines = [0 1 0 0 0 -1]'; % direction cosines for coronal positioning
        case 'sagittal'
            sag_dir_cosines = [1 0 0 0 1 0]'; % direction cosines for sagittal positioning
    end

    commandwindow();
    dir_loc = uigetdir(pwd,'Select the folder containing the sagittal localiser image'); % Select folder containing sagittal localiser
    addpath(genpath(dir_loc)); % Add to path

    Files = [dir(fullfile(dir_loc,'*.dcm'));dir(fullfile(dir_loc,'*.'));dir(fullfile(dir_loc,'*.IMA'))]; % List all files within .dcm extension
    Files = {Files.name}'; % Convert struct to cell
    Files = char(setdiff(Files,{'.','..'}));
    if isempty(Files)
        Files = ls(dir_loc);
        Files = char(setdiff(Files,{'.','..'})); % Remove upper and lower directories
    end

    if size(Files,1) > 1
        waitfor(msgbox('Multiple images detected in localiser folder. Only retrieving first sagittal image.'));
        for k = 1:size(Files,1)
            temp_obj = images.internal.dicom.DICOMFile(deblank(char(Files(k,:))));
            dir_cosines(:,k) = temp_obj.getAttributeByName('ImageOrientationPatient'); % output direction cosines
        end
        sag_ind = find(sum(dir_cosines == sag_dir_cosines)==6,1); % Find first image with matching direction cosines
        [img_loc,obj_loc,~] = ACR_OpenDICOM(Files(sag_ind,:));
    else
        [img_loc,obj_loc,~] = ACR_OpenDICOM(Files);
    end
else
    img_loc = [];
    obj_loc = [];
    dir_loc = [];
end

switch options.Orientation
    case 'coronal'
        img_loc = fliplr(imrotate(img_loc,-90,'bilinear','crop')); % reorient image
    case 'sagittal'
        img_loc = fliplr(imrotate(img_loc,-90,'bilinear','crop')); % reorient image
end
%% Series ACR 
answer = questdlg('How many ACR series would you like to open?', ...
	'Axial Datasets', ...
	'One','Two','Two');

switch answer
    case 'One'
        answer = 1; % Convert text to number
    case 'Two'
        answer = 2;
end

for m = 1:answer
    commandwindow();
    dir_ACR = uigetdir(pwd,'Select the folder containing the ACR series'); % Select folder containing ACR series
    addpath(genpath(dir_ACR)); % Add to path
    Files = [dir(fullfile(dir_ACR,'*.dcm'));dir(fullfile(dir_ACR,'*.'));dir(fullfile(dir_ACR,'*.IMA'))]; % List all files within .dcm extension

    Files = {Files.name}'; % Convert struct to cell
    Files = char(setdiff(Files,{'.','..'}));
    if isempty(Files)
        Files = dir(fullfile(dir_ACR));
        Files = {Files.name}';
        Files = char(setdiff(Files,{'.','..'})); % Remove upper and lower directories
    end

    z = zeros(1,size(Files,3)); % Initialise slice location array
 
    if size(Files,1) == 1 % Multi-frame check
        img_ACR(:,:,:,m) = squeeze(double(dicomread(deblank(char(Files(m,:))))));
        obj_ACR = images.internal.dicom.DICOMFile(deblank(char(Files(m,:))));
        list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');

        z = [list.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient(3)...
            list.Item_2.PlanePositionSequence.Item_1.ImagePositionPatient(3)...
            list.Item_3.PlanePositionSequence.Item_1.ImagePositionPatient(3)...
            list.Item_4.PlanePositionSequence.Item_1.ImagePositionPatient(3)...
            list.Item_5.PlanePositionSequence.Item_1.ImagePositionPatient(3)...
            list.Item_6.PlanePositionSequence.Item_1.ImagePositionPatient(3)...
            list.Item_7.PlanePositionSequence.Item_1.ImagePositionPatient(3)...
            list.Item_8.PlanePositionSequence.Item_1.ImagePositionPatient(3)...
            list.Item_9.PlanePositionSequence.Item_1.ImagePositionPatient(3)...
            list.Item_10.PlanePositionSequence.Item_1.ImagePositionPatient(3)...
            list.Item_11.PlanePositionSequence.Item_1.ImagePositionPatient(3)];

        [~,ind] = sort(z); % Sort slice locations into sequential order

        img_ACR(:,:,:,m) = img_ACR(:,:,ind,m); % Reorder based on ascending slice location

        orientation = round(list.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient);

        if orientation == [0 1 0 0 0 -1]' % SAGITTAL
            img_ACR = imrotate(img_ACR,-90,'bilinear','crop');
            zero_check_1 = length(find(0.1*max(img_ACR(:,:,1),[],'all') < img_ACR(:,:,1)==0));
            zero_check_11 = length(find(0.1*max(img_ACR(:,:,11),[],'all') < img_ACR(:,:,11)==0));
            if zero_check_1 < zero_check_11
                img_ACR = flipdim(img_ACR,3);
            end
        end
    else
        [img_ACR(:,:,:,m),obj_ACR(m),z] = ACR_OpenDICOM(Files);
        
        [~,ind] = sort(z); % Sort slice locations into sequential order

        img_ACR(:,:,:,m) = img_ACR(:,:,ind,m); % Reorder based on ascending slice location

        check = ACR_SliceInversionCheck(img_ACR(:,:,:,m),obj_ACR(m));
        if check == 1
            img_ACR(:,:,:,m) = flipdim(img_ACR(:,:,:,m),3);
        end
        try
            switch obj_ACR.getAttributeByName('0051,100E')
                case 'Sag'
                    img_ACR(:,:,:,m) = imrotate(img_ACR(:,:,:,m),-90,'bilinear','crop');
            end
        catch
            orientation = round(obj_ACR(m).getAttributeByName('ImageOrientationPatient'),3);

            if orientation == [0 1 0 0 0 -1]' % SAGITTAL
                img_ACR(:,:,:,m) = imrotate(img_ACR(:,:,:,m),-90,'bilinear','crop');
            end
        end
    end
end




