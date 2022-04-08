% by Yassine Azma (Mar 2022)
%
% This script takes the DICOM metadata object and determines whether the
% pre-scan normalisation filter is enabled.

% function query = ACR_PreScanNormCheck(obj_ACR)

ManufacturerStr = obj_ACR.getAttributeByName('Manufacturer');

if contains(ManufacturerStr,'Philips','IgnoreCase',true)
    Manufacturer = 'Philips';
elseif contains(ManufacturerStr,'Siemens','IgnoreCase',true)
    Manufacturer = 'Siemens';
elseif contains(ManufacturerStr,'GE','IgnoreCase',true)
    Manufacturer = 'GE';
end

switch Manufacturer
    case 'Siemens'
        if isempty(obj_ACR.getAttributeByName('ImageType')) % Multi-frame check
            list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');
            query = contains(list.Item_1.Private_2005_140f.Item_1.ImageType,'Norm','IgnoreCase',true);
        else
            ImageType = obj_ACR.getAttributeByName('ImageType');
            query = contains(ImageType,'Norm','IgnoreCase',true);
        end
    case 'Philips'
        if isempty(obj_ACR.getAttributeByName('RescaleType')) % Multi-frame check
            list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');
            query = contains(list.Item_1.Private_2005_140f.Item_1.RescaleType,'normalized');
        else
            query = contains(obj_ACR.getAttributeByName('RescaleType'),'normalized');
        end
end
    