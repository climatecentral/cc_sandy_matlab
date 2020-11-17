%Run this script after running "WarpSandy" to build the
%sandyWarpedModelCompleteData matrix, which contains latitude and longitude
%in the first two columns, and the corrected max water heights for
%simulations alt01-alt09 (alt06 being historical) in the remaining columns
    warpIdx = 6;
    errorFunc = reshape(sandyWarp(6,:,:),resY,resX);
    
    numSims = 9;
    numModelPoints = size(sandyModelCompleteData,1);
    sandyWarpedModelCompleteData = sandyModelCompleteData;
    
    rad = 0.02 * warpIdx;
    
    warpInterpolant = scatteredInterpolant(x(~isnan(errorFunc)), y(~isnan(errorFunc)), errorFunc(~isnan(errorFunc)));
    distInterpolant = griddedInterpolant(x', y', distances');
    
    modelPointsDist = distInterpolant(sandyModelCompleteData(:,1), sandyModelCompleteData(:,2));
    modelWarp = warpInterpolant(sandyModelCompleteData(:,1), sandyModelCompleteData(:,2));

    modelPointsToWarp = find(modelPointsDist < rad);
    
    modelBelowWater = sandyModelCompleteData;
    warpedModelBelowWater = sandyModelCompleteData;
    for sim = 1:numSims
        completeDataCol = sim+3;
        sandyWarpedModelCompleteData(modelPointsToWarp,completeDataCol) = sandyModelCompleteData(modelPointsToWarp,completeDataCol) - modelWarp(modelPointsToWarp);

    end