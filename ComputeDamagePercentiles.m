function [ damagePercentiles, deltaGMSL, naturalSL, historicalSL ] = ComputeDamagePercentiles( sandyTable, naturalCalibrationStruct, historicalCalibrationStruct, percentiles, stagger )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
    
            
    naturalSLMatrix = naturalCalibrationStruct.sl;
    naturalYears = naturalCalibrationStruct.time;

    historicalSLMatrix = historicalCalibrationStruct.sl;
    historicalYears = historicalCalibrationStruct.time;

    yearOfInterest=2012;
    naturalColsIndicesOfInterest = naturalYears == yearOfInterest;
    historicalColsIndicesOfInterest = historicalYears == yearOfInterest;

    naturalSL = naturalSLMatrix(:,naturalColsIndicesOfInterest);
    historicalSL = historicalSLMatrix(:,historicalColsIndicesOfInterest);


    deltaGMSL = (naturalSL - historicalSL) ./ 100;


    damages = SandyDamagesGivenDeltaGMSL(sandyTable, deltaGMSL);
    damagePercentiles = prctile(damages, 100-percentiles, 2);

end

