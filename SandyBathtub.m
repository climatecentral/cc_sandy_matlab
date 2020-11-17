stationSandyHeight = [1.795, 1.915, 3.438, 1.691, 2.836, 1.876];

stationLat = [38.97, 39.356, 40.7, 41.05, 41.17, 41.36];
stationLon = [-74.96, -74.42, -74.01, -71.96, -73.18, -72.09];


adjustments = [-0.1,-0.06,-0.02,0.02,0.06,0.10,-0.14,-0.04,0.00] - 0.10;

bathtubInterpolant = scatteredInterpolant(stationLon', stationLat', stationSandyHeight', 'nearest','nearest');

    
    count = 0;

    try
        parpool(9);
    catch
        
    end

    for tileLatIdx = 1:length(tileLatList)
        tileLat = tileLatList(tileLatIdx);

        for tileLonIdx = 1:length(tileLonList)
            count = count + 1;

            disp([num2str(count) '/' num2str(length(tileLatList) * length(tileLonList))]);

            tileLon = tileLonList(tileLonIdx);
            
            formatSpec = '%10.1f\n';
            R = georasterref('Latlim',[tileLat, tileLat + tileWidth],'Lonlim',[tileLon, tileLon + tileWidth],'RasterSize',[tileRes tileRes],'ColumnsStartFrom','north');

            tileFile = ['n' num2str(tileLat,formatSpec) 'w' num2str(-tileLon,formatSpec) '.tif'];


            [tileX, tileY] = meshgrid(linspace(tileLon,tileLon+tileWidth,tileRes/10),linspace(tileLat+tileWidth,tileLat,tileRes/10));

            pointDist = distInterpolant(tileX, tileY);

        usablePoints = pointDist < rad;

        if sum(usablePoints) == 0

            continue
        end
        
             bathtubGeotiff = bathtubInterpolant(tileX, tileY);
             bathtubGeotiff = resizem(bathtubGeotiff, [tileRes, tileRes], 'nearest');
            parfor sim=1:9 
                
                outFile = ['/data/slr1/ss2/lidar/sandy/geotiffs/bathtub/alt0' num2str(sim) '/' tileFile];
                
                geotiffwrite(outFile, single(bathtubGeotiff + adjustments(sim)), R);
            end
        end
        
    end
    