function [ output_args ] = PreprocessCounterfactualData_blocks( )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    
    upperStates = {'CT', 'NJ', 'NY'};
    

    simOrder = [6, 5, 4, 9, 3, 8, 2, 1, 7];
    
    tableDeltaGMSL = {'x_0_00m', 'x_0_04m', 'x_0_08m', 'x_0_10m', 'x_0_12m', 'x_0_14m', 'x_0_16m', 'x_0_20m ','x_0_24m'};
    
    
    for stateNum = 1:length(upperStates)
       state = upperStates{stateNum};
       
       disp(state);
       
       blockTable = readtable(['/data/slr1/ss2/lidar/blocks/' state '.unel.11m.csv']);
        
       propDamageTable = table();
       floodHeightTable = table();
       landTable = table();
       
       propDamageTable.ID = blockTable.GEOID10;
       propDamageTable.VARIABLE = cell(height(propDamageTable), 1);
       propDamageTable.VARIABLE(:) = {'PropDamage'};
       propDamageTable{:, 3:7} = 0;
       
       for simNum = simOrder
           propDamageTable.(['alt0' num2str(simNum)]) = zeros(height(propDamageTable), 1);
       end
       
       floodHeightTable = propDamageTable;
       floodHeightTable.VARIABLE(:) = {'FloodHeight'};
       
       landTable = propDamageTable;
       landTable.VARIABLE(:) = {'Land'};
       
       populationTable = propDamageTable;
       populationTable.VARIABLE(:) = {'Land'};
       
       
       for simNum = simOrder
           
          damageFile  = readtable(['/data/slr1/ss2/lidar/blocks/' state '.hans_damage.tidelcontiglevees.sandy_alt0' num2str(simNum) '_warp'], 'format', '%f%f', 'FileType', 'text');
          floodHeightFile  = readtable(['/data/slr1/ss2/lidar/blocks/' state '.averagedepth.land.sandy_alt0' num2str(simNum) '_warp'], 'format', '%f%f', 'FileType', 'text');
   
          landFile = readtable(['/data/slr1/ss2/lidar/blocks/' state '.area.land.tidelcontiglevees.sandy_alt0' num2str(simNum) '_warp'], 'format', '%f%f', 'FileType', 'text');
      
          
          for rowNum = 1:height(landFile)

            break;
            blockID = landFile.Var1(rowNum);
            landVal = landFile.Var2(rowNum);

            if blockID < 0
                continue
            end
            landTable.land(blockID+1) = landVal;
          end

          
          for rowNum = 1:height(damageFile)
             blockID = damageFile.Var1(rowNum);
             damageVal = damageFile.Var2(rowNum);
             
             
             rowInMainTable = find(blockTable.GEOID10 == blockID);
             
             propDamageTable.(['alt0' num2str(simNum)])(rowInMainTable) = damageVal;
          end
          for rowNum = 1:height(floodHeightFile)
             blockID = floodHeightFile.Var1(rowNum);
             floodHeightVal = floodHeightFile.Var2(rowNum);
             
             
             rowInMainTable = find(blockTable.FID == blockID);
             
             floodHeightTable.(['alt0' num2str(simNum)])(rowInMainTable) = floodHeightVal;
          end
    
    
       end

       rowsToDelete = (propDamageTable.alt06(:)==0 | floodHeightTable.alt06(:) == 0);
       propDamageTable(rowsToDelete, :) = [];
       floodHeightTable(rowsToDelete, :) = [];
       
       disp('Done Loading Blocks');
    
       sandyTable = [propDamageTable;floodHeightTable];
       
       alt06 = sandyTable.alt06(:);
       
        damagesTable = sandyTable(:, 1:7);
        deltaGMSLTable = table();

        deltaGMSLTable.Model = cell(9,1);
        deltaGMSLTable.Calibration = cell(9,1);
       
        LoadNaturalCMIP35Index;
        LoadHistoricalCMIP35Index;


        deltaGMSLTable = table();

        deltaGMSLTable.Model = cell(9,1);
        deltaGMSLTable.Calibration = cell(9,1);


        counterfactualDir = '/fs/data/bengis/Papers/SandyWithoutSLR/Bitterman counterfactuals/';


        models = {'CMIP5', 'LinRate', 'Mean'};
        calibrations = {'Mn', 'Mar', 'Mean'};

        counterfactualStruct = struct();

        for modelNum = 1:length(models)
           for calibrationNum = 1:length(calibrations)
               counterfactualStruct.(models{modelNum}).(calibrations{calibrationNum}).deltaGMSL = [];
               counterfactualStruct.(models{modelNum}).(calibrations{calibrationNum}).damagePercentiles = [];
           end
        end
    
        cmip5Num = 0;

        percentiles = [5 17 50 83 95];

        for percentileNum = 1:length(percentiles)
           damagesTable.(['p' num2str(percentiles(percentileNum))]) = zeros(height(damagesTable), 1); 
           deltaGMSLTable.(['p' num2str(percentiles(percentileNum))]) = zeros(9, 1);
        end
        
        damagesTable.alt06 = alt06;  
        

        %Heuristic
        for modelNum=2:3
           model = models{modelNum};



           for calibrationNum=1:2
               calibration = calibrations{calibrationNum};

               disp([model ' / ' calibration]);

               naturalFile = [counterfactualDir '/' calibration '_' model '.mat'];

               load(naturalFile);
               naturalStruct = S;

               historicalFile = [counterfactualDir '/' calibration '_HadCRUT.mat'];
               load(historicalFile);
               historicalStruct = S;

               [ damagePercentiles, deltaGMSL ] = ComputeDamagePercentiles(sandyTable, naturalStruct, historicalStruct, percentiles, 1); 

               counterfactualStruct.(model).(calibration).deltaGMSL = deltaGMSL;
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

            for calibrationNum = 1:2

                calibration = calibrations{calibrationNum};


                naturalCalibrationStruct = naturalStruct.(calibration);
                historicalCalibrationStruct = historicalStruct.(calibration);

                [ damagePercentiles, deltaGMSL ] = ComputeDamagePercentiles(sandyTable, naturalCalibrationStruct, historicalCalibrationStruct, percentiles, 1); 

                counterfactualStruct.CMIP5.(calibrations{calibrationNum}).deltaGMSL(:,cmip5Num) = deltaGMSL;
                counterfactualStruct.CMIP5.(calibrations{calibrationNum}).damagePercentiles(:,:,cmip5Num) = damagePercentiles;

            end
        end

        for calibrationNum = 1:2
            counterfactualStruct.CMIP5.(calibrations{calibrationNum}).deltaGMSL = round(mean(counterfactualStruct.CMIP5.(calibrations{calibrationNum}).deltaGMSL, 2),3);
            counterfactualStruct.CMIP5.(calibrations{calibrationNum}).damagePercentiles = round(mean(counterfactualStruct.CMIP5.(calibrations{calibrationNum}).damagePercentiles, 3), 3);
        end

        for modelNum=1:3
           model=models{modelNum};
           counterfactualStruct.(model).(calibrations{3}).deltaGMSL = (counterfactualStruct.(model).Mn.deltaGMSL + counterfactualStruct.(model).Mar.deltaGMSL) ./ 2;
           counterfactualStruct.(model).(calibrations{3}).damagePercentiles = (counterfactualStruct.(model).Mn.damagePercentiles + counterfactualStruct.(model).Mar.damagePercentiles) ./ 2;

           for calibrationNum=1:3
               calibration = calibrations{calibrationNum};

               deltaGMSLTableRow = (modelNum - 1) * 3 + calibrationNum;

               deltaGMSLTable.Model{deltaGMSLTableRow} = model;
               deltaGMSLTable.Calibration{deltaGMSLTableRow} = calibration;
               deltaGMSLTable{deltaGMSLTableRow, 3:end} = prctile(counterfactualStruct.(model).(calibration).deltaGMSL, percentiles);

               for percentileNum = 1:length(percentiles)
                  percentile = percentiles(percentileNum);

                  percentileStr = ['p' num2str(percentile)];
                  deltaPercentileStr = ['delta_' percentileStr];
                  fracPercentileStr = ['frac_' percentileStr];

                  damagesTable.(percentileStr)(:) = counterfactualStruct.(model).(calibration).damagePercentiles(:,percentileNum);
                  damagesTable.(deltaPercentileStr) = damagesTable.alt06(:) - damagesTable.(percentileStr)(:);
                  
                  damagesTable.(deltaPercentileStr)(damagesTable.(deltaPercentileStr)(:) < 0) = 0;
                  damagesTable.(fracPercentileStr) = damagesTable.(deltaPercentileStr)(:) ./ damagesTable.alt06(:);
                  
               end
               damageTableCondensed = table();
               
               propDamageSubTable = damagesTable(strcmp(damagesTable.VARIABLE, 'PropDamage'),{'ID','frac_p50'});
               floodHeightSubTable = damagesTable(strcmp(damagesTable.VARIABLE, 'FloodHeight'),{'ID','p50'});
               
               damageTableCondensed.ID = propDamageSubTable.ID;
               damageTableCondensed.FloodHeight_p50 = floodHeightSubTable.p50(:);
               damageTableCondensed.PropDamage_frac_p50 = propDamageSubTable.frac_p50(:);
               
               damageTableCondensed(damageTableCondensed.PropDamage_frac_p50 == 0 | damageTableCondensed.FloodHeight_p50 == 0,:) = [];
               
               writetable(damagesTable, ['output/' state '_sandy_cfdamages_' model '_' calibration '.blocks.csv']);
               writetable(damageTableCondensed, ['output/' state '_sandy_cfdamages_' model '_' calibration '.blocks.condensed.csv']);
           end
        end
    end
end

