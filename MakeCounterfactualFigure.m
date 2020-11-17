function [ output_args ] = MakeCounterfactualFigure( )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C, 'EdgeColor','none', 'FaceAlpha', 0.25);

    
    run('LoadNaturalCMIP35Index.m');
    run('LoadHistoricalCMIP35Index.m');
    sandyTable = readtable('/data/slr1/ss2/lidar/results/sandy_v2.csv');
    
    alt06 = sandyTable.('x_0_00m')(:);
    
    damagesTable = sandyTable(:, 1:7);
    deltaGMSLTable = table();
    
    deltaGMSLTable.Model = cell(9,1);
    deltaGMSLTable.Calibration = cell(9,1);
    
    
    counterfactualDir = '/fs/data/bengis/Papers/SandyWithoutSLR/Bitterman counterfactuals/';
    
    
    models = {'CMIP5',  'LinRate', 'Mean'};
    %models = {'CMIP5'};
    calibrations = {'Mn', 'Mar', 'Mean'};
        
    counterfactualStruct = struct();
    
    obsTable = readtable([counterfactualDir '/GMSL_Hay.csv']);
    
    obsSLR = [obsTable.GMSL(:)' - obsTable.std90(:)'; obsTable.GMSL(:)'; obsTable.GMSL(:)' + obsTable.std90(:)']; 
    
    for modelNum = 1:length(models)
       for calibrationNum = 1:length(calibrations)
           counterfactualStruct.(models{modelNum}).(calibrations{calibrationNum}).deltaGMSL = [];
       end
    end
    
    cmip5Num = 0;
    
    percentiles = [5 50 95];
    
    for percentileNum = 1:length(percentiles)
       damagesTable.(['p' num2str(percentiles(percentileNum))]) = zeros(height(damagesTable), 1); 
       deltaGMSLTable.(['p' num2str(percentiles(percentileNum))]) = zeros(9, 1);
       
       damagesTable.(['delta_p' num2str(percentiles(percentileNum))]) = zeros(height(damagesTable), 1);
    end
    
    damagesTable.alt06 = alt06;  
    
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
           
           
     
           counterfactualStruct.(model).(calibration).naturalGMSL = (naturalStruct.sl);
           counterfactualStruct.(model).(calibration).historicalGMSL = (historicalStruct.sl);
           counterfactualStruct.(model).(calibration).deltaGMSL = (historicalStruct.sl - naturalStruct.sl);
            
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
            
            counterfactualStruct.CMIP5.(calibration).naturalGMSL(:,:,cmip5Num) = (naturalCalibrationStruct.sl);
            counterfactualStruct.CMIP5.(calibration).historicalGMSL(:,:,cmip5Num) = (historicalCalibrationStruct.sl);
            
            counterfactualStruct.CMIP5.(calibration).deltaGMSL(:,:,cmip5Num) = (historicalCalibrationStruct.sl - naturalCalibrationStruct.sl);
            
        end
        
    end
    
    
    for modelNum = 1:length(models)
        model = models{modelNum};
        
        %counterfactualStruct.(model).Mean.naturalGMSL = (counterfactualStruct.(model).Mn.naturalGMSL + counterfactualStruct.(model).Mar.naturalGMSL) ./ 2;
        %counterfactualStruct.(model).Mean.historicalGMSL = (counterfactualStruct.(model).Mn.historicalGMSL + counterfactualStruct.(model).Mar.historicalGMSL) ./ 2;
        %counterfactualStruct.(model).Mean.deltaGMSL = (counterfactualStruct.(model).Mn.deltaGMSL + counterfactualStruct.(model).Mar.deltaGMSL) ./ 2;
        
        counterfactualStruct.(model).Mean.naturalGMSL = cat(3,counterfactualStruct.(model).Mn.naturalGMSL, counterfactualStruct.(model).Mar.naturalGMSL);
        counterfactualStruct.(model).Mean.historicalGMSL = cat(3,counterfactualStruct.(model).Mn.historicalGMSL, counterfactualStruct.(model).Mar.historicalGMSL);
        counterfactualStruct.(model).Mean.deltaGMSL = cat(3,counterfactualStruct.(model).Mn.deltaGMSL, counterfactualStruct.(model).Mar.deltaGMSL);
        
        
        %Take overall means
        for calibrationNum = 1:3

            calibration = calibrations{calibrationNum};

            counterfactualStruct.(model).(calibration).naturalGMSL(:,114:end,:) = [];
            counterfactualStruct.(model).(calibration).historicalGMSL(:,114:end,:) = [];
            counterfactualStruct.(model).(calibration).deltaGMSL(:,114:end,:) = [];

            counterfactualStruct.(model).(calibration).naturalGMSL = mean(counterfactualStruct.(model).(calibration).naturalGMSL, 3);
            counterfactualStruct.(model).(calibration).historicalGMSL = mean(counterfactualStruct.(model).(calibration).historicalGMSL, 3);
            counterfactualStruct.(model).(calibration).deltaGMSL = mean(counterfactualStruct.(model).(calibration).deltaGMSL, 3);
            
            counterfactualStruct.(model).(calibration).naturalGMSL = prctile(counterfactualStruct.(model).(calibration).naturalGMSL, percentiles, 1);
            counterfactualStruct.(model).(calibration).historicalGMSL = prctile(counterfactualStruct.(model).(calibration).historicalGMSL, percentiles, 1);
            counterfactualStruct.(model).(calibration).deltaGMSL = prctile(counterfactualStruct.(model).(calibration).deltaGMSL, percentiles, 1);
        end


        
    end
    years = 1900:2012;
    fig = figure('Position', [0 0 1200 700], 'Color', [1 1 1]);
    boundingBox = [1900 2010 -4 23];
    
    
    
    
    subplot(1,3,1); hold on; box on;
    
    
    %plot(years, obsSLR(1, :), 'k--');
    %plot(years, obsSLR(3, :), 'k--');
    %fill_between_lines(years, obsSLR(1, :), obsSLR(3, :), [0 0 0]);
    
    
    orangeH = plot(years, counterfactualStruct.Mean.Mean.historicalGMSL(2, :), '-', 'Color', [1 0.5 0], 'LineWidth', 2);
    %plot(years, counterfactualStruct.Mean.Mean.historicalGMSL(1, :), 'b--');
    %plot(years, counterfactualStruct.Mean.Mean.historicalGMSL(3, :), 'b--');
    fill_between_lines(years, counterfactualStruct.Mean.Mean.historicalGMSL(1, :), counterfactualStruct.Mean.Mean.historicalGMSL(3, :), [1 0.5 0]);
    
    magentaH = plot(years, counterfactualStruct.CMIP5.Mean.historicalGMSL(2, :), 'm', 'LineWidth', 2);
    %plot(years, counterfactualStruct.CMIP5.Mean.historicalGMSL(1, :), 'r--');
    %plot(years, counterfactualStruct.CMIP5.Mean.historicalGMSL(3, :), 'r--');
    fill_between_lines(years, counterfactualStruct.CMIP5.Mean.historicalGMSL(1, :), counterfactualStruct.CMIP5.Mean.historicalGMSL(3, :), [1 0 1]);
    
    
    blackH = plot(1900:2010, obsSLR(2, :), 'k', 'LineWidth', 2);
    
    axis(boundingBox);
    title('Observed/Historical SLR'); xlabel('Year'); ylabel('SLR (cm)');
    
    subplot(1,3,2); hold on; box on;
    
    blueH = plot(years, counterfactualStruct.Mean.Mean.naturalGMSL(2, :), 'b-', 'LineWidth', 2);
    %plot(years, counterfactualStruct.Mean.Mean.naturalGMSL(1, :), 'b--');
    %plot(years, counterfactualStruct.Mean.Mean.naturalGMSL(3, :), 'b--');
    fill_between_lines(years, counterfactualStruct.Mean.Mean.naturalGMSL(1, :), counterfactualStruct.Mean.Mean.naturalGMSL(3, :), [0 0 1]);
    
    
    greenH = plot(years, counterfactualStruct.LinRate.Mean.naturalGMSL(2, :), '-', 'Color', [0 0.5 0], 'LineWidth', 2);
    %plot(years, counterfactualStruct.LinRate.Mean.naturalGMSL(1, :), '--', 'Color', [0 0 0.5], 'LineWidth', 2);
    %plot(years, counterfactualStruct.LinRate.Mean.naturalGMSL(3, :), '--', 'Color', [0 0 0.5], 'LineWidth', 2);
    fill_between_lines(years, counterfactualStruct.LinRate.Mean.naturalGMSL(1, :), counterfactualStruct.LinRate.Mean.naturalGMSL(3, :), [0 0.5 0]);
    
    
    plot(years, counterfactualStruct.CMIP5.Mean.naturalGMSL(2, :), 'm', 'LineWidth', 2);
    %plot(years, counterfactualStruct.CMIP5.Mean.naturalGMSL(1, :), 'r--');
    %plot(years, counterfactualStruct.CMIP5.Mean.naturalGMSL(3, :), 'r--');
    fill_between_lines(years, counterfactualStruct.CMIP5.Mean.naturalGMSL(1, :), counterfactualStruct.CMIP5.Mean.naturalGMSL(3, :), [1 0 1]);
    axis(boundingBox);
    title('Counterfactual'); xlabel('Year'); ylabel('SLR (cm)');
    
    subplot(1,3,3); hold on; box on;
    
    plot(years, counterfactualStruct.Mean.Mean.deltaGMSL(2, :), 'b', 'LineWidth', 2);
    %plot(years, counterfactualStruct.Mean.Mean.deltaGMSL(1, :), 'b--');
    %plot(years, counterfactualStruct.Mean.Mean.deltaGMSL(3, :), 'b--');
    fill_between_lines(years, counterfactualStruct.Mean.Mean.deltaGMSL(1, :), counterfactualStruct.Mean.Mean.deltaGMSL(3, :), [0 0 1]);
    
    plot(years, counterfactualStruct.LinRate.Mean.deltaGMSL(2, :), '-', 'Color', [0 0.5 0], 'LineWidth', 2);
    %plot(years, counterfactualStruct.LinRate.Mean.deltaGMSL(1, :), '--', 'Color', [0 0 0.5], 'LineWidth', 2);    
    %plot(years, counterfactualStruct.LinRate.Mean.deltaGMSL(3, :), '--', 'Color', [0 0 0.5], 'LineWidth', 2);
    fill_between_lines(years, counterfactualStruct.LinRate.Mean.deltaGMSL(1, :), counterfactualStruct.LinRate.Mean.deltaGMSL(3, :), [0 0.5 0]);
    
    
    plot(years, counterfactualStruct.CMIP5.Mean.deltaGMSL(2, :), 'm', 'LineWidth', 2);
    %plot(years, counterfactualStruct.CMIP5.Mean.deltaGMSL(1, :), 'r--');
    %plot(years, counterfactualStruct.CMIP5.Mean.deltaGMSL(3, :), 'r--');
    fill_between_lines(years, counterfactualStruct.CMIP5.Mean.deltaGMSL(1, :), counterfactualStruct.CMIP5.Mean.deltaGMSL(3, :), [1 0 1]);
    axis(boundingBox);
    title('Difference'); xlabel('Year'); ylabel('SLR (cm)');
    
    legend([blackH, magentaH, orangeH, blueH, greenH], {'Observed', 'CMIP5', 'HadCRUT4',  'Stable', 'Cooling'}, 'Orientation', 'horizontal');
    
    save('cf_data.mat','obsSLR','counterfactualStruct');
end

