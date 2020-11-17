datadir = 'E:/data/slr1/ss2/lidar/sandy/geotiffs/';

outFileNoWarp = @(simNum,filename)[datadir '/july_simulations/alt0' num2str(simNum) '/no_warp/' filename];
outFileWarp = @(simNum,filename)[datadir '/july_simulations/alt0' num2str(simNum) '/warp/' filename];

tileRes = 2222;

tileWidth = 0.1;
tileMinLat = floor(minLat * 10) / 10;
tileMinLon = floor(minLon * 10) / 10;
tileMaxLat = floor(maxLat * 10) / 10;
tileMaxLon = floor(maxLon * 10) / 10;

numSims = 9;

tileLatList = tileMinLat:tileWidth:tileMaxLat;
tileLonList = tileMinLon:tileWidth:tileMaxLon;

if ~exist('modelFloodInterpolant')
    modelFloodInterpolant = cell(numSims * 2, 1);
    
    for sim=1:numSims
        disp(sim);
        mkdir([datadir '/july_simulations/alt0' num2str(sim) '/no_warp']);
        mkdir([datadir '/july_simulations/alt0' num2str(sim) '/warp']);

        modelFloodInterpolant{sim * 2 - 1} = scatteredInterpolant(sandyModelCompleteData(:,1),sandyModelCompleteData(:,2), sandyModelCompleteData(:,sim+3), 'nearest','nearest');
        modelFloodInterpolant{sim * 2} = scatteredInterpolant(sandyWarpedModelCompleteData(:,1),sandyWarpedModelCompleteData(:,2), sandyWarpedModelCompleteData(:,sim+3), 'nearest','nearest');

    end
end

try
   parpool(numSims*2); 
catch

end

parfor_progress(length(tileLatList) * length(tileLonList));

count = 0;

for tileLatIdx = 1:length(tileLatList)
    tileLat = tileLatList(tileLatIdx);

    for tileLonIdx = 1:length(tileLonList)
        count = count + 1;
        
        disp([num2str(count) '/' num2str(length(tileLatList) * length(tileLonList))]);

        tileLon = tileLonList(tileLonIdx);
        
        formatSpec = '%10.1f\n';
        R = georasterref('Latlim',[tileLat, tileLat + tileWidth],'Lonlim',[tileLon, tileLon + tileWidth],'RasterSize',[tileRes tileRes],'ColumnsStartFrom','north');

        tileFile = ['n' num2str(tileLat,formatSpec) 'w' num2str(-tileLon,formatSpec) '.tif'];
        

        [tileX, tileY] = meshgrid(linspace(tileLon,tileLon+tileWidth,tileRes),linspace(tileLat+tileWidth,tileLat,tileRes));
        pointDist = distInterpolant(tileX, tileY);

        usablePoints = pointDist < rad;

        if sum(usablePoints) == 0

            continue
        end

        tileX = tileX(usablePoints);
        tileY = tileY(usablePoints);
        parfor run=1:numSims*2
            
            tileData = zeros(tileRes, tileRes);
            
            sim = ceil(run / 2);
            
            filenameWarp = outFileWarp(sim,tileFile);
            filenameNoWarp = outFileNoWarp(sim,tileFile);
            tileData(usablePoints) = modelFloodInterpolant{run}(tileX,tileY);
            
            tileData(isnan(tileData)) = -9999;
            if mod(run, 2) == 0
                geotiffwrite(filenameWarp,single(tileData),R);
            else
                geotiffwrite(filenameNoWarp,single(tileData),R);
            end
        end

    end
end
