baselineSim = 'alt06';

sandyModelCompleteData_Scaled = sandyModelCompleteData;
sandyWarpedModelCompleteData_Scaled = sandyWarpedModelCompleteData;

adjustments = [-0.1,-0.06,-0.02,0.02,0.06,0.10,-0.14,-0.04,0.00] - 0.10;

for sim=1:9
    sandyModelCompleteData_Scaled(:,3+sim) = sandyModelCompleteData(:,3+6) + adjustments(sim);
    sandyWarpedModelCompleteData_Scaled(:,3+sim) = sandyWarpedModelCompleteData_Scaled(:,3+6) + adjustments(sim);
    
    warpTypes = {'warp','no_warp'};
    
    disp(sim);
    disp(['Bias: ' num2str(mean(sandyModelCompleteData_Scaled(:,3+sim) - sandyModelCompleteData(:,3+sim)))]);
    disp(['RMSE: ' num2str(sqrt(mean((sandyModelCompleteData_Scaled(:,3+sim) - sandyModelCompleteData(:,3+sim)).^2)))]);

    for warpType = warpTypes
        baseDir = ['/data/slr1/ss2/lidar/sandy/geotiffs/july_simulations/alt06/' warpType{:}];
        scaledBaseDir = ['/data/slr1/ss2/lidar/sandy/geotiffs/scaled/alt0' num2str(sim) '/' warpType{:}];
        
        try
            mkdir(scaledBaseDir);
        catch
        end
        
        fileList = dir([baseDir '/*.tif']);
        numFiles = numel(fileList);
        
        for fileNum = 1:numFiles
            
            file = fileList(fileNum).name;
            
            
            inFullFile = [baseDir '/' file];
            outFullFile = [scaledBaseDir '/' file];
            
            disp([num2str(fileNum) '/' num2str(numFiles) '(alt0' num2str(sim) ' ' warpType{:} ')']);
    
            
            [tile, R] = geotiffread(inFullFile);
            
            tile = tile + adjustments(sim);
            tile(tile < -9999) = -9999;
            
            geotiffwrite(outFullFile, tile, R);
            
        end
    end
end
