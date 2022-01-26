function z = ACR_RetrieveSliceThickness(obj_ACR)

if isempty(obj_ACR.getAttributeByName('SliceThickness')) % Multi-frame check
    list = obj_ACR.getAttributeByName('PerFrameFunctionalGroupsSequence');
    z = list.Item_1.PixelMeasuresSequence.Item_1.SliceThickness;
else
    z = obj_ACR.getAttributeByName('SliceThickness'); % retrieve ACR prescribed slice thickness
end