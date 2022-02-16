function [img,obj,z] = ACR_OpenDICOM(Files)

for k = 1:size(Files,1)
    temp = double(dicomread(deblank(char(Files(k,:)))));
    obj = images.internal.dicom.DICOMFile(deblank(char(Files(k,:))));
    Intercept = obj.getAttributeByName('RescaleIntercept');
    Slope = obj.getAttributeByName('RescaleSlope');
    if size(Files,1) > 1 
        z(k) =  obj.getAttributeByName('SliceLocation'); % Extract slice location
    else
        z = [];
    end
    % Rescaling
    if isempty(Intercept) || isempty(Slope) % If either are empty, do not attempt rescaling
        img(:,:,k) = temp;
    else
        if strcmpi(obj.getAttributeByName('Manufacturer'),'Philips')
            SS = obj.getAttributeByName('Private_2005_100E');

            if isempty(SS) || sum(size(SS))>0
                img(:,:,k) = temp;
            else
                img(:,:,k) = (temp + (Intercept/Slope))./SS;
            end
        else
            img(:,:,k) = Slope.*temp + Intercept;
        end
    end
end