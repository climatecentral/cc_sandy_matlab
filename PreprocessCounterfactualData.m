function [ output_args ] = PreprocessCounterfactualData( )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    LoadNaturalCMIP35Index;
    LoadHistoricalCMIP35Index;
    
    slCItoStr = @(slPercentiles) [ num2str(slPercentiles(1)) ' - ' num2str(slPercentiles(3)) ];
               
    
    sandyTable = readtable('/data/slr1/ss2/lidar/results/sandy_v2.csv');
    
    alt06 = sandyTable.('x_0_00m')(:);
    
    damagesTable = sandyTable(:, 1:7);
    summaryTable = table();
    
    summaryTable.Model = cell(9,1);
    summaryTable.Calibration = cell(9,1);
    
    counterfactualDir = '/fs/data/bengis/Papers/SandyWithoutSLR/Bitterman counterfactuals/';
    
    
    models = {'CMIP5', 'LinRate', 'Mean'};
    %models = {'CMIP5'};
    calibrations = {'Mn', 'Mar', 'Mean'};
        
    cmip5Table = table();
    
    counterfactualStruct = struct();
    
    for modelNum = 1:length(models)
       for calibrationNum = 1:length(calibrations)
           counterfactualStruct.(models{modelNum}).(calibrations{calibrationNum}).deltaGMSL = [];
           counterfactualStruct.(models{modelNum}).(calibrations{calibrationNum}).damagePercentiles = [];
       end
    end
    
    cmip5Num = 0;
    
    percentiles = [5 50 95];
    
    for percentileNum = 1:length(percentiles)
       damagesTable.(['p' num2str(percentiles(percentileNum))]) = zeros(height(damagesTable), 1); 
       
       damagesTable.(['delta_p' num2str(percentiles(percentileNum))]) = zeros(height(damagesTable), 1);
       damagesTable.(['frac_p' num2str(percentiles(percentileNum))]) = zeros(height(damagesTable), 1);
    end
    
    summaryTable.deltaSL_median = zeros(9, 1);
    summaryTable.deltaSL_CI = cell(9, 1);
    summaryTable.historicalSL_median = zeros(9, 1);
    summaryTable.historicalSL_CI = cell(9, 1);
    summaryTable.naturalSL_median = zeros(9, 1);
    summaryTable.naturalSL_CI = cell(9, 1);
    
    
    damagesTable.alt06 = alt06;  
    
    %Heuristic
    for modelNum=2:3
       
       try
           model = models{modelNum};
       catch
           continue
       end
       
        
       for calibrationNum=1:2
           calibration = calibrations{calibrationNum};
           
           disp([model ' / ' calibration]);
           
           naturalFile = [counterfactualDir '/' calibration '_' model '.mat'];
            
           load(naturalFile);
           naturalStruct = S;
        
           historicalFile = [counterfactualDir '/' calibration '_HadCRUT.mat'];
           load(historicalFile);
           historicalStruct = S;
           
           [ damagePercentiles, deltaGMSL, naturalSL, historicalSL ] = ComputeDamagePercentiles(sandyTable, naturalStruct, historicalStruct, percentiles, 1); 
            
           counterfactualStruct.(model).(calibration).deltaGMSL = deltaGMSL;
           counterfactualStruct.(model).(calibration).naturalGMSL = naturalSL;
           counterfactualStruct.(model).(calibration).historicalGMSL = historicalSL;
           counterfactualStruct.(model).(calibration).damagePercentiles = damagePercentiles;
            
       end
       
    end
    
    %CMIP5
    for naturalRowNum = 1:height(naturalCMIP35index)
        
        disp([num2str(naturalRowNum) '/' num2str(height(naturalCMIP35index))]);
        naturalIdx = naturalCMIP35index.Ind(naturalRowNum);
        
        model = naturalCMIP35index.exofmodelswhos{naturalRowNum};
        
        if strcmp(model, 'GFDL-CM3')
            continue
        end
                
        replicant = naturalCMIP35index.edatais{naturalRowNum};
        
        matchingHistoricalRow = find(strcmp(historicalCMIP35index.exofmodelswhos(:), model) & strcmp(historicalCMIP35index.edatais(:), replicant));
        
        if isempty(matchingHistoricalRow)
            continue
        end
        
        
        historicalIdx = historicalCMIP35index.Ind(matchingHistoricalRow);
        
        naturalFile = [counterfactualDir '/CMIPs_extended_thru_2012/CMIP5_hist_&_nat50yrTrend/CF-CMIP5_' num2str(naturalIdx) '.mat'];
        historicalFile = [counterfactualDir '/CMIPs_extended_thru_2012/CMIP5_hist_&_nat50yrTrend/CMIP5-hist_' num2str(historicalIdx) '.mat'];
        
        load(naturalFile);
        naturalStruct = S;
        
        load(historicalFile);
        historicalStruct = S;
        
        cmip5Num = cmip5Num + 1;
        
        counterfactualStruct.CMIP5.names{cmip5Num} = [model ', ' replicant];
        
        
        for calibrationNum = 1:2
           
            calibration = calibrations{calibrationNum};
            
            
            naturalCalibrationStruct = naturalStruct.(calibration);
            historicalCalibrationStruct = historicalStruct.(calibration);
            
            [ damagePercentiles, deltaGMSL, naturalSL, historicalSL ] = ComputeDamagePercentiles(sandyTable, naturalCalibrationStruct, historicalCalibrationStruct, percentiles, 1); 
            
            counterfactualStruct.CMIP5.(calibrations{calibrationNum}).historicalGMSL_individual(:,cmip5Num) = historicalSL;
            counterfactualStruct.CMIP5.(calibrations{calibrationNum}).naturalGMSL_individual(:,cmip5Num) = naturalSL;
            
            
            counterfactualStruct.CMIP5.(calibrations{calibrationNum}).deltaGMSL_individual(:,cmip5Num) = deltaGMSL;
            counterfactualStruct.CMIP5.(calibrations{calibrationNum}).damagePercentiles_individual(:,:,cmip5Num) = damagePercentiles;
            
        end
        
    end
    
    for calibrationNum = 1:2
        counterfactualStruct.CMIP5.(calibrations{calibrationNum}).historicalGMSL = round(mean(counterfactualStruct.CMIP5.(calibrations{calibrationNum}).historicalGMSL_individual, 2),3);
        counterfactualStruct.CMIP5.(calibrations{calibrationNum}).naturalGMSL = round(mean(counterfactualStruct.CMIP5.(calibrations{calibrationNum}).naturalGMSL_individual, 2),3);
        counterfactualStruct.CMIP5.(calibrations{calibrationNum}).deltaGMSL = round(mean(counterfactualStruct.CMIP5.(calibrations{calibrationNum}).deltaGMSL_individual, 2),3);
        counterfactualStruct.CMIP5.(calibrations{calibrationNum}).damagePercentiles = round(mean(counterfactualStruct.CMIP5.(calibrations{calibrationNum}).damagePercentiles_individual,3),2);
    end
    
    for modelNum=1:length(models)
       model=models{modelNum};
       
       if strcmp(model, 'CMIP5')
           counterfactualStruct.(model).Mean.naturalGMSL_individual = (counterfactualStruct.(model).Mn.naturalGMSL_individual + counterfactualStruct.(model).Mar.naturalGMSL_individual) ./ 2;
           counterfactualStruct.(model).Mean.historicalGMSL_individual = (counterfactualStruct.(model).Mn.historicalGMSL_individual + counterfactualStruct.(model).Mar.historicalGMSL_individual) ./ 2;        
           counterfactualStruct.(model).Mean.deltaGMSL_individual = (counterfactualStruct.(model).Mn.deltaGMSL_individual + counterfactualStruct.(model).Mar.deltaGMSL_individual) ./ 2;
           counterfactualStruct.(model).Mean.damagePercentiles_individual = (counterfactualStruct.(model).Mn.damagePercentiles_individual + counterfactualStruct.(model).Mar.damagePercentiles_individual) ./ 2;
       
       end
       
       counterfactualStruct.(model).Mean.historicalGMSL = (counterfactualStruct.(model).Mn.historicalGMSL + counterfactualStruct.(model).Mar.historicalGMSL) ./ 2;
       counterfactualStruct.(model).Mean.naturalGMSL = (counterfactualStruct.(model).Mn.naturalGMSL + counterfactualStruct.(model).Mar.naturalGMSL) ./ 2;
       counterfactualStruct.(model).Mean.deltaGMSL = (counterfactualStruct.(model).Mn.deltaGMSL + counterfactualStruct.(model).Mar.deltaGMSL) ./ 2;
       counterfactualStruct.(model).Mean.damagePercentiles = (counterfactualStruct.(model).Mn.damagePercentiles + counterfactualStruct.(model).Mar.damagePercentiles) ./ 2;
       
       for calibrationNum=1:3
           calibration = calibrations{calibrationNum};
           
           summaryTableRow = (modelNum - 1) * 3 + calibrationNum;
           
           historicalSL = round(squeeze(prctile(counterfactualStruct.(model).(calibration).historicalGMSL, [5 50 95], 1)'),1);
           naturalSL = round(squeeze(prctile(counterfactualStruct.(model).(calibration).naturalGMSL, [5 50 95], 1)'),1);
           deltaSL = round(squeeze(prctile(-counterfactualStruct.(model).(calibration).deltaGMSL*100, [5 50 95], 1)'),1);
           
           
           summaryTable.Model{summaryTableRow} = model;
           summaryTable.Calibration{summaryTableRow} = calibration;
           summaryTable.naturalSL_median(summaryTableRow) = naturalSL(2);
           summaryTable.naturalSL_CI{summaryTableRow} = slCItoStr(naturalSL);
           
           summaryTable.historicalSL_median(summaryTableRow) = historicalSL(2);
           summaryTable.historicalSL_CI{summaryTableRow} = slCItoStr(historicalSL);
           
           summaryTable.deltaSL_median(summaryTableRow) = deltaSL(2);
           summaryTable.deltaSL_CI{summaryTableRow} = slCItoStr(deltaSL);
           
          
           for percentileNum = 1:length(percentiles)
              percentile = percentiles(percentileNum);
              
              percentileStr = ['p' num2str(percentile)];
              deltaPercentileStr = ['delta_' percentileStr];
              fracPercentileStr = ['frac_' percentileStr];
              
              damagesTable.(percentileStr)(:) = counterfactualStruct.(model).(calibration).damagePercentiles(:,percentileNum);
              damagesTable.(deltaPercentileStr)(:) = damagesTable.alt06(:) - damagesTable.(percentileStr)(:);
              
              damagesTable.(deltaPercentileStr)(damagesTable.(deltaPercentileStr)(:) < 0) = 0;
              damagesTable.(fracPercentileStr)(:) = damagesTable.(deltaPercentileStr)(:) ./ damagesTable.alt06(:);
           end
           
           if strcmp(model, 'CMIP5')
               
               cmip5Table.Name = counterfactualStruct.CMIP5.names';
               
               for calibrationNum = 1:length(calibrations)
                    calibration = calibrations{calibrationNum};
                  
                    historicalSL = round(squeeze(prctile(counterfactualStruct.CMIP5.(calibration).historicalGMSL_individual, [5 50 95], 1)'),1);
                    naturalSL = round(squeeze(prctile(counterfactualStruct.CMIP5.(calibration).naturalGMSL_individual, [5 50 95], 1)'),1);
                    deltaSL = round(squeeze(prctile(-100*counterfactualStruct.CMIP5.(calibration).deltaGMSL_individual, [5 50 95], 1)'),1);
                    
                    cmip5Table.([calibration '_historical_median']) = historicalSL(:,2);
                    cmip5Table.([calibration '_cf_median']) = naturalSL(:,2);
                    cmip5Table.([calibration '_diff_median']) = deltaSL(:,2);
                    
                    for simNum = 1:size(historicalSL, 1)
                        cmip5Table.([calibration '_historical_CI']){simNum} = slCItoStr(historicalSL(simNum,:));
                        cmip5Table.([calibration '_cf_CI']){simNum} = slCItoStr(naturalSL(simNum,:));
                        cmip5Table.([calibration '_diff_CI']){simNum} = slCItoStr(deltaSL(simNum,:));
                    end
                    
               end
               writetable(cmip5Table, ['output/sandy_indiv_cmips.csv']);
           end
           writetable(damagesTable, ['output/sandy_cfdamages_' model '_' calibration '.csv']);
       end
    end
    
    writetable(summaryTable, 'output/sandy_summary.csv');
end

