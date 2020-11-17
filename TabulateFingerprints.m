function [outputArg1,outputArg2] = TabulateFingerprints(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    inDir = 'E:/data/slr1/ss2/lidar/sandy/fingerprints';
    
    gridFiles = FilesInDirectory([inDir '/*.grd']);
    
    stationNames = {'Halifax', 'Eastport', 'Bar Harbor', 'Portland', 'Nantucket Island', 'Woods Hole', 'Newport', 'New London', 'Montauk', 'New York (The Battery)', 'Sandy Hook', 'Annapolis (naval Academy)', 'Solomon''s Island (Biol. lab.)', 'Swells point, Hampton Roads'};
    
    stationLons = [-63.58,-66.98,-68.21,-70.25, -70.10,-70.68,-71.33, -72.09, -71.96, -74.01, -74.00, -76.48, -76.45, -76.33];
    stationLats = [44.67, 44.90, 44.39, 43.66, 41.29, 41.52, 41.51, 41.36, 41.05, 40.70, 40.46, 38.98, 38.32, 36.85];
    
    nycLon = -74.006;
    nycLat = 40.7128;
    
    maskNCID = netcdf.open([inDir '/land_mask/slm_720.grd']);
    maskGrid = netcdf.getVar(maskNCID, 2);
    
    maskLons = netcdf.getVar(maskNCID, 0);
    maskLats = netcdf.getVar(maskNCID,1);
    maskGrid = flipud(maskGrid');
    
    colsToFlip = maskLons > 180;
    
    maskGrid = [maskGrid(:,colsToFlip) maskGrid(:,~colsToFlip)];
    
    pixelDegreeWidths = 0.5;
    
    pixelLonDistKmByLat = arrayfun(@(x) lldistkm([maskLats(x),0], [maskLats(x),pixelDegreeWidths]), 1:length(maskLats));
    
    maskWeights = maskGrid .* pixelLonDistKmByLat';
    
    for fileNum = 5:length(gridFiles)
        
        fileName = gridFiles{fileNum};
        disp(fileName);
        
        mapLabel = strrep(fileName, '.grd', '');
        ncid = netcdf.open([inDir '/' fileName], 'NC_NOWRITE');
        gridYears = netcdf.getVar(ncid, 0);
        
        rslID = netcdf.inqVarID(ncid,'RSL');
        disp(rslID);
        if ~exist('outTable', 'var')
            outTable = table();
            outTable.Year = gridYears;
        end
        
        gridLats = netcdf.getVar(ncid, 1);
        gridLons = netcdf.getVar(ncid, 2);
        
        if max(gridLats)>90
            temp = gridLats;
            gridLats = gridLons;
            gridLons = temp;
        end
        
        gridLats = flipud(gridLats);
       % gridLons = gridLons + 180;
        
        gridLons = gridLons + 180;
        pixToFlip = gridLons>360;
        
        gridLons = [gridLons(pixToFlip)-540; gridLons(~pixToFlip)-180];
        
        gridData = netcdf.getVar(ncid, rslID);
        
        gridData = permute(gridData, [2 1 3]);
        gridData = flip(gridData, 1);
        gridData = cat(2, gridData(:,pixToFlip,:), gridData(:,~pixToFlip,:));
        
        [meshGridLons, meshGridLats] = meshgrid(gridLons, gridLats);
        
        
        nycTimeSeries = zeros(length(gridYears), 1);
        globalAvgTimeSeries = zeros(length(gridYears), 1);
        nycFingerprintCoeff = zeros(length(gridYears), 1);
        
        for yearNum = 1:length(gridYears)
           
            gridDataForYear = gridData(:,:,yearNum);
            
            globalAvg = sum(sum(gridDataForYear.*maskWeights)) / sum(maskWeights(:));
            nycVal = interp2(meshGridLons, meshGridLats,  gridDataForYear, single(nycLon), single(nycLat));
        
            nycTimeSeries(yearNum) = nycVal;
            globalAvgTimeSeries(yearNum) = globalAvg;
            nycFingerprintCoeff(yearNum) = nycVal / globalAvg;
        end
        
        outTable.([mapLabel '_nyc_grd_val']) = nycTimeSeries;
        outTable.([mapLabel '_global_avg_grd_val']) = globalAvgTimeSeries;
        outTable.([mapLabel '_nyc_fingerprint_coeff']) = nycFingerprintCoeff;
        
    end
    writetable(outTable, 'output/fingerprints_data/nyc_fingerprints.csv');
end

