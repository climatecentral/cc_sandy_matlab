if ~exist('blocksShp')
    blocksShp.NJ = shaperead('E:/data/slr1/ss2/lidar/blocks/NJ.unel.11m.shp', 'UseGeoCoords', true);
    blocksShp.NY = shaperead('E:/data/slr1/ss2/lidar/blocks/NY.unel.11m.shp', 'UseGeoCoords', true);
    blocksShp.CT = shaperead('E:/data/slr1/ss2/lidar/blocks/CT.unel.11m.shp', 'UseGeoCoords', true);
end

lonlim = [-74.4 ,-73.4];
latlim = [40.2,41.2];


states = {'NJ', 'NY', 'CT'};

doDamageCalc = 1;


for stateNum=1:length(states)
    state = states{stateNum};
    disp(state)
    damageTable = readtable(['output/' state '_sandy_cfdamages_CMIP5_Mean.csv']);
    
    damageModeled = damageTable{strcmp(damageTable.VARIABLE, 'PropDamage'), [1 10]};
    damageHistorical = damageTable{strcmp(damageTable.VARIABLE, 'PropDamage'), [1 end]};
    
    depthModeled = damageTable{strcmp(damageTable.VARIABLE, 'FloodHeight'), [1 10]};
    depthHistorical = damageTable{strcmp(damageTable.VARIABLE, 'FloodHeight'), [1 end]};
    
    if doDamageCalc
        blocksToUse.(state) = [];
        
        for i=1:length(blocksShp.(state))
           blocksShp.(state)(i).damageHistorical = 0;
           blocksShp.(state)(i).damageModeled = 0;

           blocksShp.(state)(i).depthHistorical = 0;
           blocksShp.(state)(i).depthModeled = 0;
           blocksShp.(state)(i).depthDiff = 0;
           geoid = str2num(blocksShp.(state)(i).GEOID10);

           blockLat = blocksShp.(state)(i).Lat;
           blockLon = blocksShp.(state)(i).Lon;

           if ~(any(blockLat>latlim(1)) && any(blockLat < latlim(2)) && ...
                any(blockLon>lonlim(1)) && any(blockLon < lonlim(2)))
            continue
           end

           damageHistoricalIdx = find(damageHistorical(:,1) == geoid);
           damageModeledIdx = find(damageModeled(:,1) == geoid);


           depthHistoricalIdx = find(depthHistorical(:,1) == geoid);
           depthModeledIdx = find(depthModeled(:,1) == geoid);

           if ~isempty(depthHistoricalIdx)
               blocksShp.(state)(i).depthHistorical = depthHistorical(depthHistoricalIdx, 2);   
           end
           if ~isempty(depthModeledIdx)
               blocksShp.(state)(i).depthModeled = depthModeled(depthModeledIdx, 2);   
           end 

           blocksShp.(state)(i).depthDiff =  1 - (blocksShp.(state)(i).depthModeled / blocksShp.(state)(i).depthHistorical);



           if ~isempty(damageHistoricalIdx)
               blocksShp.(state)(i).damageHistorical = (damageHistorical(damageHistoricalIdx, 2));   
           end

           if ~isempty(damageModeledIdx)
               blocksShp.(state)(i).damageModeled = (damageModeled(damageModeledIdx, 2));   
           end 

           if blocksShp.(state)(i).depthHistorical == 0
               continue
           end

           blocksShp.(state)(i).damageDiff = 1 -(blocksShp.(state)(i).damageModeled/blocksShp.(state)(i).damageHistorical);

           %blocksShp.(state)(i).depthDiff = min(blocksShp.(state)(i).depthDiff, 0.3);
           %blocksShp.(state)(i).damageDiff = min(blocksShp.(state)(i).damageDiff, 0.3);
           
           %blocksShp.(state)(i).damageHistorical = min(blocksShp.(state)(i).damageHistorical, 1000);
           
           if blocksShp.(state)(i).depthHistorical > 0 && blocksShp.(state)(i).damageDiff > 0 && blocksShp.(state)(i).depthDiff > 0 && blocksShp.(state)(i).damageHistorical > 2
               blocksToUse.(state) = [blocksToUse.(state) i];
           else
               blocksShp.(state)(i).damageDiff = -1;
               blocksShp.(state)(i).depthDiff = -1;
               
           end
        end
    end
end

minDamage = 0;
maxDamage = 1;

%minDamage = 0;
%maxDamage = (1000);

damageColorMap = colormap(jet(100));

minDepth = 0;
maxDepth =1;
depthColorMap =  colormap(jet(100));

damageDiffColors = makesymbolspec('Polygon',  {'Default','EdgeColor','none','FaceColor','none'},{'damageDiff',[(minDamage) (maxDamage)],'FaceColor',damageColorMap});
depthDiffColors = makesymbolspec('Polygon',  {'Default','EdgeColor','none','FaceColor','none'},{'depthDiff',[(minDepth) (maxDepth)],'FaceColor',depthColorMap});

fig = figure('Units','inches','position', [0, 0, 20, 10]);

njCoast = shaperead('/data/slr1/ss2/lidar/states/NJ.shoreline.simple.high.shp', 'UseGeoCoords', true);
nyCoast = shaperead('/data/slr1/ss2/lidar/states/NY.shoreline.simple.high.shp', 'UseGeoCoords', true);
ctCoast = shaperead('/data/slr1/ss2/lidar/states/CT.shoreline.simple.high.shp', 'UseGeoCoords', true);

oceanColor = [0.5 0.7 0.9];

set(fig,'PaperPosition', [0 0 7 7]);


subtightplot(1,2,1);

ax = usamap(latlim,lonlim);
set(ax, 'Visible', 'off')
setm(ax, 'FFaceColor', oceanColor);
gridm off; mlabel off; plabel off


title('Property Damage Differencce (0-100%)'); 
mapshow(njCoast, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',1)
mapshow(nyCoast, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',1)
mapshow(ctCoast, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',1)

mapshow(blocksShp.NJ(blocksToUse.NJ), 'DisplayType', 'polygon', 'SymbolSpec', damageDiffColors);
mapshow(blocksShp.NY(blocksToUse.NY), 'DisplayType', 'polygon', 'SymbolSpec', damageDiffColors);
mapshow(blocksShp.CT(blocksToUse.CT), 'DisplayType', 'polygon', 'SymbolSpec', damageDiffColors);


subtightplot(1,2,2);
ax = usamap(latlim,lonlim);
set(ax, 'Visible', 'off')
setm(ax, 'FaceColor', oceanColor);
gridm off; mlabel off; plabel off

title('Anthropogenic Water Height %'); 
mapshow(njCoast, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',1)
mapshow(nyCoast, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',1)
mapshow(ctCoast, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',1)

mapshow(blocksShp.NJ(blocksToUse.NJ), 'DisplayType', 'polygon', 'SymbolSpec', depthDiffColors);
mapshow(blocksShp.NY(blocksToUse.NY), 'DisplayType', 'polygon', 'SymbolSpec', depthDiffColors);
mapshow(blocksShp.CT(blocksToUse.CT), 'DisplayType', 'polygon', 'SymbolSpec', depthDiffColors);

