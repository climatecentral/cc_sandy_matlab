
if ~exist('sandyObservationData')
    sandyStationsToRemove = readtable('data/station_to_remove.csv');
    sandyObservationData = readtable('data/2012_Sandy_AC2_YesDatumCorrectedNoGeoidOffset.csv');
    
    idsToRemove = [];
% 
    for i=1:numel(sandyObservationData.OBS_NAME)

        if strcmp(sandyObservationData.OBS_TYPE{i}, 'HWM')

            stationID=sandyObservationData.OBS_NAME{i};
        else
            stationID=sandyObservationData.OBS_ID{i};

        end
        if any(strcmp(stationID,sandyStationsToRemove.STATION_ID(:)))
           idsToRemove = [idsToRemove, i]; 
           disp(stationID)

        end
    end

    idsToRemove = unique([idsToRemove find(sandyObservationData.LON_DEG_EAST < -73.84 & sandyObservationData.LAT_DEG_NORTH > 41.0)']);
    idsToRemove = unique([idsToRemove find(not(cellfun('isempty',strfind(sandyObservationData.OBS_DESCRIPTION,'Poor'))))']);
    sandyObservationData(idsToRemove,:) = [];

end
    
if ~exist('sandyModelMesh') || ~exist('sandyModelCompleteData')
    [sandyModelMesh,~] = read_adcirc_mesh('data/simulations_july_2016/research_Sandy_01_fromNACP43_Rivers.grd/Sandy_01_fromNACP43_Rivers.grd');
    numSims = 9;
    sandyModelCompleteData = zeros(numel(sandyModelMesh.x), numSims);
    sandyModelCompleteData(:,1) = sandyModelMesh.x;
    sandyModelCompleteData(:,2) = sandyModelMesh.y;
    sandyModelCompleteData(:,3) = sandyModelMesh.z;
    for i = 1:numSims
       simFileName = ['data/simulations_july_2016/research_TP_0001_HIS_Tides_1_SLC_0_RFC_1_maxele_alt0' num2str(i) '.63/NACCS_TP_0001_HIS_Tides_1_SLC_0_RFC_1_maxele.63'];
        disp(simFileName);
       [~,~,sandyModelCompleteData(:, 3+i)] = read_63(simFileName,1);
       sandyModelCompleteData(sandyModelCompleteData(:,3+i)==-99999, 3+i)=nan;
    end
    
    modelNodesToRemove = [];
end


conversionData = csvread('data/NAVDtoMSL.csv',1,0);



if ~exist('modelNodesToRemove') || length(modelNodesToRemove) == 0
    modelNodesToRemove = sandyModelMesh.x * 0;
    minLon = min(sandyObservationData.LON_DEG_EAST) - 0.1;
    maxLon = max(sandyObservationData.LON_DEG_EAST) + 0.1;
    minLat = min(sandyObservationData.LAT_DEG_NORTH) - 0.1;
    maxLat = max(sandyObservationData.LAT_DEG_NORTH) + 0.1;

    modelNodesToRemove(sandyModelMesh.x < minLon,:) = 1;
    modelNodesToRemove(sandyModelMesh.x > maxLon,:) = 1;
    modelNodesToRemove(sandyModelMesh.y < minLat,:) = 1;
    modelNodesToRemove(sandyModelMesh.y > maxLat,:) = 1;

    sandyModelCompleteData(find(modelNodesToRemove),:) = [];
end
%sandyObservationData(idsToRemove,:) = [];

conversionLat = conversionData(:,5);
conversionLon = conversionData(:,4);

conversionFactor = conversionData(:,6);
conversionFactorInterpolant = scatteredInterpolant(conversionLon,conversionLat,conversionFactor, 'linear','nearest');

sandyErrorTable = sandyObservationData(:,:);

sandyObservationLat = sandyErrorTable.LAT_DEG_NORTH;
sandyObservationLon = sandyErrorTable.LON_DEG_EAST;
sandyObservationElev =  (sandyErrorTable.OBS_METERS_NAVD88);

modelLon = sandyModelCompleteData(:,1);
modelLat = sandyModelCompleteData(:,2);

for i=1:numSims
    disp(i);
    modelNanNodes = isnan(sandyModelCompleteData(:,3+i));
    
    modelHeightInterpolant = scatteredInterpolant(modelLon(~modelNanNodes),modelLat(~modelNanNodes),sandyModelCompleteData(~modelNanNodes, 3 + i), 'linear','nearest');
    %sandyModelCompleteData(modelNanNodes, 3+i) = -9999;%modelHeightInterpolant(modelLon(modelNanNodes),modelLat(modelNanNodes));
    
    sandyModelCompleteData(:, 3+i) = sandyModelCompleteData(:, 3+i) - conversionFactorInterpolant(modelLon, modelLat);
    
    sandyModelElevAtObservedPoints = modelHeightInterpolant(sandyObservationLon, sandyObservationLat) - conversionFactorInterpolant(sandyObservationLon, sandyObservationLat);

    sandyError =  sandyModelElevAtObservedPoints - sandyObservationElev;

    sandyErrorTable.(['MODEL_0' num2str(i) '_MAXELEV_NAVD88']) = sandyModelElevAtObservedPoints;
    sandyErrorTable.(['CC_MODEL_0' num2str(i) '_MINUS_OBS']) = sandyError;
end