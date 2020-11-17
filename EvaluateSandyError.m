function [ sandyError ] = EvaluateSandyError( observationData, modelData, modelDataCol )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    modelNaN = isnan(modelData(:,modelDataCol));
    modelHeightInterpolant = scatteredInterpolant(modelData(~modelNaN,1),modelData(~modelNaN,2),modelData(~modelNaN,modelDataCol));
    

    obsLon = observationData.LON_DEG_EAST;
    obsLat = observationData.LAT_DEG_NORTH;
    obsElev = observationData.OBS_METERS_NAVD88;
    
    sandyModelElevAtObservedPoints = modelHeightInterpolant(obsLon, obsLat);
    
    
    sandyError =  [obsLon obsLat (sandyModelElevAtObservedPoints - obsElev)];
    
    sandyError(obsElev<=0,:) = [];
    obsElev(obsElev<=0) = [];
    disp(['Bias: ' num2str(mean(sandyError(:,3)))]);
    disp(['RMSE: ' num2str(sqrt(mean(sandyError(:,3).^2))) ]);
    disp(['NRMSE: ' num2str(sqrt(mean(sandyError(:,3).^2)) / mean(obsElev)) ]);

end

