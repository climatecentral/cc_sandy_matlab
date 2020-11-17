function [outputArg1,outputArg2] = MakeDamageShapefiles()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    states = {'CT', 'NJ','NY'};
    %states = {'NJ','NY'};
    damage_alt06 = readtable('data/damage/blocks_depth.alt06.warp.simulation.csv', 'Format', '%s%f');
    damage_alt09 = readtable('data/damage/blocks_depth.alt09.warp.simulation.csv', 'Format', '%s%f');
    damage_alt05 = readtable('data/damage/blocks_depth.alt05.warp.simulation.csv', 'Format', '%s%f');
    damage_alt01 = readtable('data/damage/blocks_depth.alt01.warp.simulation.csv', 'Format', '%s%f');
    
    damage_alt06.CB_ = cellfun(@(x) x(2:end), damage_alt06.CB_, 'UniformOutput', false);
    damage_alt09.CB_ = cellfun(@(x) x(2:end), damage_alt09.CB_, 'UniformOutput', false);
    damage_alt05.CB_ = cellfun(@(x) x(2:end), damage_alt05.CB_, 'UniformOutput', false);
    damage_alt01.CB_ = cellfun(@(x) x(2:end), damage_alt01.CB_, 'UniformOutput', false);
    
    for stateNum = 1:length(states)
        state = states{stateNum};
        disp(state);
        
        stateShp = shaperead(['E:/data/slr1/ss2/lidar/sandy/data_for_hans/' state '.shp']);
        
        
        blockIDs = {stateShp.GEOID10};
        parfor_progress(floor(length(blockIDs) / 1000));
        
        
        results = zeros(length(blockIDs),14);
        parfor blockNum = 1:length(blockIDs)
            if mod(blockNum, 1000) == 0
                parfor_progress;
            end
            blockID = blockIDs{blockNum};
            
            alt06Damage = 0;
            alt09Damage = 0;
            alt05Damage = 0;
            alt01Damage = 0;
            
            
            
            try
                alt06Damage = damage_alt06.TOTAL(strcmp(damage_alt06.CB_(:), blockID));
            catch
            end
            
            try
                alt09Damage = damage_alt09.TOTAL(strcmp(damage_alt09.CB_(:), blockID));
            catch 
            end
            
            try
                alt05Damage = damage_alt05.TOTAL(strcmp(damage_alt05.CB_(:), blockID));
            catch 
            end
            
            try
                alt01Damage = damage_alt01.TOTAL(strcmp(damage_alt01.CB_(:), blockID));
            catch 
            end
            
            
            if isempty(alt06Damage)
                alt06Damage = 0; 
            end
             
            if isempty(alt09Damage)
                alt09Damage = 0;
            end
             
            if isempty(alt05Damage)
                alt05Damage = 0;
            end
             
            if isempty(alt01Damage)
                alt01Damage = 0;
            end
            
            alt06DamagePerSqKm = alt06Damage / (stateShp(blockNum).ALAND10/1000000);
            alt09DamagePerSqKm = alt09Damage / (stateShp(blockNum).ALAND10/1000000);
            alt05DamagePerSqKm = alt05Damage / (stateShp(blockNum).ALAND10/1000000);
            alt01DamagePerSqKm = alt01Damage / (stateShp(blockNum).ALAND10/1000000);
            
            if isnan(alt06DamagePerSqKm) || isinf(alt06DamagePerSqKm)
                alt06DamagePerSqKm = 0;
            end
            
            if isnan(alt09DamagePerSqKm) || isinf(alt09DamagePerSqKm)
                alt09DamagePerSqKm = 0;
            end
            if isnan(alt05DamagePerSqKm) || isinf(alt05DamagePerSqKm)
                alt05DamagePerSqKm = 0;
            end
            if isnan(alt01DamagePerSqKm) || isinf(alt01DamagePerSqKm)
                alt01DamagePerSqKm = 0;
            end
            
            damageDiffPerSqKmAlt09 = alt06DamagePerSqKm - alt09DamagePerSqKm;
            damageDiffPerSqKmAlt09Log10 = log10(damageDiffPerSqKmAlt09);
            
            
            damageDiffPerSqKmAlt05 = alt06DamagePerSqKm - alt05DamagePerSqKm;
            damageDiffPerSqKmAlt05Log10 = log10(damageDiffPerSqKmAlt05);
            
            damageDiffPerSqKmAlt01 = alt06DamagePerSqKm - alt01DamagePerSqKm;
            damageDiffPerSqKmAlt01Log10 = log10(damageDiffPerSqKmAlt01);
            
            
            if damageDiffPerSqKmAlt09Log10 < 0 || ~isreal(damageDiffPerSqKmAlt09Log10) || isnan(damageDiffPerSqKmAlt09Log10) || isinf(damageDiffPerSqKmAlt09Log10)
                damageDiffPerSqKmAlt09Log10 = 0;
            end
            
            if damageDiffPerSqKmAlt05Log10 < 0 || ~isreal(damageDiffPerSqKmAlt05Log10) || isnan(damageDiffPerSqKmAlt05Log10) || isinf(damageDiffPerSqKmAlt05Log10)
                damageDiffPerSqKmAlt05Log10 = 0;
            end
            
            if damageDiffPerSqKmAlt01Log10 < 0 || ~isreal(damageDiffPerSqKmAlt01Log10) || isnan(damageDiffPerSqKmAlt01Log10) || isinf(damageDiffPerSqKmAlt01Log10)
                damageDiffPerSqKmAlt01Log10 = 0;
            end
            
            
            
            results(blockNum, :) = [alt06Damage, alt09Damage, alt05Damage, alt01Damage, alt06DamagePerSqKm, alt09DamagePerSqKm, alt05DamagePerSqKm, alt01DamagePerSqKm, damageDiffPerSqKmAlt09, damageDiffPerSqKmAlt05, damageDiffPerSqKmAlt01, damageDiffPerSqKmAlt09Log10, damageDiffPerSqKmAlt05Log10, damageDiffPerSqKmAlt01Log10 ];
            
            
        end
        
        results = num2cell(results);
        
        [stateShp.totaldam06] = results{:,1};
        [stateShp.totaldam09]=results{:,2};
        [stateShp.totaldam05]=results{:,3};
        [stateShp.totaldam01]=results{:,4};
        
        [stateShp.perkm2dam06]=results{:,5};
        [stateShp.perkm2dam09]=results{:,6};
        [stateShp.perkm2dam05]=results{:,7};
        [stateShp.perkm2dam01]=results{:,8};
        
        [stateShp.diffdam9]=results{:,9};
        [stateShp.diffdam5]=results{:,10};
        [stateShp.diffdam1]=results{:,11};
        
        [stateShp.diffdam9l10]=results{:,12};
        [stateShp.diffdam5l10]=results{:,13};
        [stateShp.diffdam1l10]=results{:,14};
        
        shapewrite(stateShp, ['data/damage/' state '.damage.shp']);
    end
end
