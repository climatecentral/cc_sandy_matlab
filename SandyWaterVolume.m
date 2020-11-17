states = {'NJ','NY','CT'};

layers = {};
layersStr = {};

for sim = 1:9
        simStr = ['alt0' num2str(sim)];


        layers{end+1} = ['bathtub/' simStr];
        
        layersStr{end+1} = ['sandy_' simStr '_bathtub'];
        
        for warp = {'warp','no_warp'}
            warp = warp{1};
            for type = {'july_simulations','scaled'}
                type = type{1};
                layers{end+1} = [type '/' simStr '/' warp];
                
                if strcmp(type, 'july_simulations')
                    layersStr{end+1} = ['sandy_' simStr '_' warp];
                else
                    layersStr{end+1} = ['sandy_' simStr '_' warp '_scaled'];
                end
            end
        end
        
end


try
   % parpool(12);
catch
end

parfor stateIdx =1:numel(states)
    
    stateStr = states{stateIdx};
    
    tileList = dir(['/data/slr0/ss2/lidar/tidelthresh/' stateStr '/sandy_alt06_bathtub/*.tif']);
    numTiles = numel(tileList);
    
    blockAverageDepthData = zeros(numel(layers),50000);
    blockInundatedCountData = zeros(numel(layers),50000);
    blockVolumeData = zeros(numel(layers),50000);
    
    for layerIdx=1:numel(layers)
        layerStr = layersStr{layerIdx};
    end
    
    for tileIdx = 1:numTiles
        tic;
       
        
        tileFile = tileList(tileIdx).name;
        [landTile, landR] = geotiffread(['/data/slr0/ss2/lidar/dryland/' stateStr '/' tileFile]);
        
        
        if ~any(landTile>0)
            continue
        end
        
        landTile = single(landTile);
        
        [elevTile, elevR] = geotiffread(['/data/slr1/ss2/lidar/unel/' stateStr '/' tileFile]);
        [blocksTile, blocksR] = geotiffread(['/data/slr0/ss2/lidar/blocks/' stateStr '/' tileFile]);
        
        blocksInTile = unique(blocksTile);
        
        lonPixelDist = distance(landR.LatitudeLimits(1), landR.LongitudeLimits(1), landR.LatitudeLimits(1), landR.LongitudeLimits(2), [6378.1,0.00335]) * 1000 / 2222; 
        latPixelDist = distance(landR.LatitudeLimits(1), landR.LongitudeLimits(1), landR.LatitudeLimits(2), landR.LongitudeLimits(1), [6378.1,0.00335]) * 1000 / 2222; 

        pixelArea = lonPixelDist * latPixelDist;
        
        if size(blockVolumeData,2) < max(blocksInTile)
            blockAverageDepthData(numel(layers),max(blocksInTile)) = 0;
            blockInundatedCountData(numel(layers),max(blocksInTile)) = 0;
            blockVolumeData(numel(layers),max(blocksInTile)) = 0;
        end
        
        for layerIdx=1:numel(layers)
            layer = layers{layerIdx};
            layerStr = layersStr{layerIdx};
            
            try
                [floodTile, floodR] = geotiffread(['/data/slr0/ss2/lidar/tidelcontiglevees/' stateStr '/' layerStr '/' tileFile]);
            catch
                continue
            end
            
            if ~any(floodTile>0)
                continue
            end
            
            pixelsOfInterest = find(floodTile == 1 & blocksTile > 0 & landTile>0)';
            
            [modelHeightTile, modelHeightR] = geotiffread(['/data/slr1/ss2/lidar/sandy/geotiffs/' layer '/' tileFile]);
            
            depthTile = (modelHeightTile - elevTile) .* landTile .* single(floodTile);
            
            %depthTile(depthTile<0) = 0;
            
            for pixelIdx=pixelsOfInterest
                
                pixelDepth = depthTile(pixelIdx);
                pixelBlock = blocksTile(pixelIdx);
                pixelWaterVolume = pixelDepth * pixelArea;
                
                
                blockVolumeData(layerIdx,pixelBlock) = blockVolumeData(layerIdx,pixelBlock) + pixelWaterVolume;
                
                prevAverage = blockAverageDepthData(layerIdx,pixelBlock);
                prevCount = blockInundatedCountData(layerIdx,pixelBlock);
                newCount = prevCount + 1;
                blockInundatedCountData(layerIdx,pixelBlock) = newCount;
                
                blockAverageDepthData(layerIdx,pixelBlock) = ((prevAverage * prevCount) + pixelDepth) / newCount;               
                
            end
        end
        disp([stateStr ' ' num2str(tileIdx) '/' num2str(numTiles) ' - ' num2str(toc)])
    end
    
    for layerIdx=1:numel(layers)
        layer = layers{layerIdx};
        layerStr = layersStr{layerIdx};
       
        outFileVolume = ['/data/slr1/ss2/lidar/blocks/' stateStr '.floodvolume.land.' layerStr];
        outFileDepth = ['/data/slr1/ss2/lidar/blocks/' stateStr '.averagedepth.land.' layerStr];
        
        outDataVolume = [(1:size(blockVolumeData,2))' blockVolumeData(layerIdx,:)'];
        outDataVolume(outDataVolume(:,2) == 0,:) = [];
        
        outDataDepth = [(1:size(blockAverageDepthData,2))' blockAverageDepthData(layerIdx,:)'];
        outDataDepth(outDataDepth(:,2) == 0,:) = [];
        
        dlmwrite(outFileVolume, outDataVolume,'precision',10);
        dlmwrite(outFileDepth, outDataDepth,'precision',10);
    end
    
end