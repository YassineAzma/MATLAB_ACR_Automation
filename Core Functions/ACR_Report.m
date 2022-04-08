function ACR_Report(L,resolvable,dz,dL,PIU,PSG,SNR,eff_res,obj_ACR,obj_loc,options)

% FOLDER FIELDS
Location = obj_ACR.getAttributeByName('InstitutionName');
Scanner = [obj_ACR.getAttributeByName('Manufacturer') ' ' obj_ACR.getAttributeByName('ManufacturerModelName') ' ',...
    num2str(obj_ACR.getAttributeByName('MagneticFieldStrength')) 'T'];

TargetFolder = [options.OutputPath '\' Location];
if not(isfolder(TargetFolder))
    mkdir(TargetFolder);
end

% SPREADSHEET FIELDS
if strcmpi(options.IncludeLocaliser,'yes')
    Date_loc = obj_loc.getAttributeByName('AcquisitionDate');
    Coil_loc = obj_loc.getAttributeByName('ReceiveCoilName');
    if isempty(Coil_loc)
        try 
            Coil_loc = obj_loc.getAttributeByName('SharedFunctionalGroupsSequence').Item_1.MRReceiveCoilSequence.Item_1.ReceiveCoilName;
        catch
            try 
                Coil_loc = obj_loc.getAttributeByName('Private_0051_100F');
            catch
                Coil_loc = [];
            end
        end
    end
    SeriesName_loc = [num2str(obj_loc.getAttributeByName('SeriesNumber')) ' - ' obj_loc.getAttributeByName('SeriesDescription')];
end

Date_ACR = obj_ACR.getAttributeByName('AcquisitionDate');
if isempty(Date_ACR)
    Date_ACR = obj_ACR.getAttributeByName('StudyDate');
end

Coil_ACR = obj_ACR.getAttributeByName('ReceiveCoilName');
if isempty(Coil_ACR)
    try
        Coil_ACR = obj_ACR.getAttributeByName('SharedFunctionalGroupsSequence').Item_1.MRReceiveCoilSequence.Item_1.ReceiveCoilName;
    catch
        try
            Coil_ACR = obj_ACR.getAttributeByName('Private_0051_100F');
        catch
            Coil_ACR = [];
        end
    end
end
% if strcmpi(ManufacturerName,'Siemens')
%     Elements = obj_ACR.getAttributeByName('Private_0021_114F');
% end
SeriesName_ACR = [num2str(obj_ACR.getAttributeByName('SeriesNumber')) ' - ' obj_ACR.getAttributeByName('SeriesDescription')];

% PROCESS DATA
if strcmpi(options.IncludeLocaliser,'yes')
    geo_acc_loc = L(1);
    GeometricAccuracyLoc = geo_acc_loc;
else
    geo_acc_loc = [];
end

geo_acc_ACR = round(mean(L(2:end)),2);
geo_acc_cov_ACR = round(100*std(L(2:end))/mean(L(2:end)),2);

if strcmpi(options.IncludeLocaliser,'yes')
    SeriesNameLoc = {SeriesName_loc};
end
SeriesNameACR = {SeriesName_ACR};

Date_ACR = [Date_ACR(end-1:end) '/' Date_ACR(end-3:end-2) '/' Date_ACR(1:4)];
Time_ACR = obj_ACR.getAttributeByName('AcquisitionTime');
if isempty(Time_ACR)
    Time_ACR = obj_ACR.getAttributeByName('StudyTime');
end
Time_ACR = [Time_ACR(1:2) ':' Time_ACR(3:4) ':' Time_ACR(5:6)];

DateTime_ACR = [Date_ACR ' - ' Time_ACR];

if strcmpi(options.IncludeLocaliser,'yes')
    Date_loc = [Date_loc(end-1:end) '/' Date_loc(end-3:end-2) '/' Date_loc(1:4)];
    Time_loc = obj_loc.getAttributeByName('AcquisitionTime');
    if isempty(Time_ACR)
        Time_loc = obj_loc.getAttributeByName('AcquisitionTime');
    end
    Time_loc = [Time_loc(1:2) ':' Time_loc(3:4) ':' Time_loc(5:6)];

    DateTime_loc = [Date_loc ' - ' Time_loc];
    Coil_loc = {Coil_loc};
end

Coil_ACR = {Coil_ACR};

GeometricAccuracyACR = geo_acc_ACR;
GeometricAccuracyCoV = geo_acc_cov_ACR;
SpatialResolution = {resolvable};
SliceThickness = round(dz,2);

if strcmpi(options.SliceDisplacementDisplay,'yes')
    SlicePosition = [num2str(round(dL(1),2)) '*, ' num2str(round(dL(2),2)) '*'];
else
    SlicePosition = [num2str(round(dL(1),2)) ', ' num2str(round(dL(2),2))];
end
Uniformity = round(PIU,2);
Ghosting = round(PSG,2);
if strcmpi(options.RicianNoiseCorrection,'yes')
    SignaltoNoise = [num2str(round(SNR,2)) '*'];
else
    SignaltoNoise = round(SNR,2);
end
MTF50 = round(eff_res,2);

Header = {'Date','Coil','Series Name','Geometric Accuracy (mm)','Coefficient of Variation (%)','Spatial Resolution'...
    ,'Slice Thickness (mm)','Slice Position (mm)','Integral Uniformity (%)','Percent Signal Ghosting (%)','SNR',...
    'Effective Resolution (mm)'};
ACR_Fields = {[DateTime_ACR,Coil_ACR,SeriesNameACR,GeometricAccuracyACR,GeometricAccuracyCoV,SpatialResolution,SliceThickness ... 
    SlicePosition,Uniformity,Ghosting,SignaltoNoise,MTF50]};

% EXISTING FIELDS HANDLING
if isfile([TargetFolder '\'  Scanner '.xls'])
    ExistTable = readtable([TargetFolder '\' Scanner '.xls'],'VariableNamingRule','Preserve');
    ExistName = table2cell(ExistTable(:,3));
    ExistDate = table2cell(ExistTable(:,1));
    CompNameACR = strcmpi(ExistName,SeriesName_ACR);
    CompDateACR = strcmpi(ExistDate,DateTime_ACR);
end
if strcmpi(options.IncludeLocaliser,'yes')
    loc_Fields = {[DateTime_loc,Coil_loc,SeriesNameLoc,GeometricAccuracyLoc cell(1,7)]};
    if isfile([TargetFolder '\' Scanner '.xls'])
        CompNameloc = strcmpi(ExistName,SeriesName_loc);
        CompDateloc = strcmpi(ExistDate,DateTime_loc);

        if sum(CompNameACR.*CompDateACR)==0 && sum(CompNameloc.*CompDateloc)==0
            T = table([loc_Fields{1,1};ACR_Fields{1,1}]);
        elseif sum(CompNameloc.*CompDateloc)==0
            T = table(loc_Fields{1,1});
        elseif sum(CompNameACR.*CompDateACR)==0
            T = table(ACR_Fields{1,1});
        else
            T = [];
        end
    else
        T = table([Header;loc_Fields{1,1};ACR_Fields{1,1}]);
    end
else
    if isfile([TargetFolder '\' Scanner '.xls'])
        if sum(CompNameACR.*CompDateACR)==0
            T = table([ACR_Fields{1,1}]);
        end
    else
        T = table([Header;ACR_Fields{1,1}]);
    end
end

if ~exist('T','var')
    T = [];
end

% WRITE TABLE
if ~isempty(T)
    writetable(T,[TargetFolder '\' Scanner '.xls'],'WriteVariableNames',0,'WriteMode','append')
    % SORT DATA
    writetable(cell2table([Header; table2cell(sortrows(readtable([TargetFolder '\' Scanner '.xls'],'VariableNamingRule','Preserve')))]),[TargetFolder '\' Scanner '.xls'],...
        'WriteVariableNames',0');
end

