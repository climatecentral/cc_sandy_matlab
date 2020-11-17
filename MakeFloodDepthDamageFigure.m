cfDamageTable = readtable('output/sandy_cfdamages_CMIP5_Mean.csv');





anthroDamageFrac = 1 - (cfDamageTable.p50(:) ./ sandyTable.x_0_00m(:));


states = {'NJ', 'NY', 'CT'};

propDamageRows = strcmp(sandyTable.VARIABLE, 'PropDamage');
floodVolumeRows = strcmp(sandyTable.VARIABLE, 'FloodVolume');
landAreaRows = strcmp(sandyTable.VARIABLE, 'Land');

correctModelRows = strcmp(sandyTable.LAYER_TYPE, 'Simulation') & strcmp(sandyTable.WARPED, 'Yes');

for stateNum = 1:length(states)
    state = states{stateNum};
   
    stateTownShp = shaperead(['/data/slr1/ss2/lidar/munis/' state '.unel.11m.shp']);
    
    for featureNum = 1:length(stateTownShp)
        ssid = [state '_Town_' stateTownShp(featureNum).GEOID10];
        
        featureRows = (strcmp(cfDamageTable.PLACE_ID, ssid));
        
        featurePropDamageFrac = anthroDamageFrac(featureRows & propDamageRows & correctModelRows);
        featureFloodVolumeFrac = anthroDamageFrac(featureRows & floodVolumeRows & correctModelRows);
        
        featureFloodHeightHistorical = sandyTable.x_0_00m(featureRows & floodVolumeRows & correctModelRows) / sandyTable.x_0_00m(featureRows & landAreaRows & correctModelRows);
        featureFloodHeightModeled = cfDamageTable.p50(featureRows & floodVolumeRows & correctModelRows) / cfDamageTable.p50(featureRows & landAreaRows & correctModelRows);
        
        featureFloodHeightAntroFrac = 1 - (featureFloodHeightModeled / featureFloodHeightHistorical);
        
       
        
    end
end