function [outputArg1,outputArg2] = MakeCountyDamageUncertaintyFigure()
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

    fontSize = 14;
    countyDamageTable = readtable('output/sandy_v2.csv');
    countyDamageTable(~strcmp(countyDamageTable.PLACE_TYPE(:), 'County') | ~strcmp(countyDamageTable.VARIABLE(:), 'PropDamage') | ~strcmp(countyDamageTable.LAYER_TYPE, 'Simulation') | ~strcmp(countyDamageTable.WARPED, 'Yes'),:) = [];
    
    
    countyDamageTable(11,:) = [];
    
    colors = zeros(height(countyDamageTable), 3);
    
    colors(strcmp(countyDamageTable.STATE(:), 'NJ'), :) = repmat([1 0 0], nnz(strcmp(countyDamageTable.STATE(:), 'NJ')), 1);
    colors(strcmp(countyDamageTable.STATE(:), 'NY'), :) = repmat([1 0 1], nnz(strcmp(countyDamageTable.STATE(:), 'NY')), 1);
    colors(strcmp(countyDamageTable.STATE(:), 'CT'), :) = repmat([0 0 1], nnz(strcmp(countyDamageTable.STATE(:), 'CT')), 1);
    
    
    medianASLR=9.6;
    lowASLR=5.6;
    highASLR=15.6;
    fig = figure('Position', [0 0 800 800]);
    hold on;
    
    countyDamageASLRs = [0 4 8 10 12 14 16 20 24];
    countyDamages = countyDamageTable{:, 8:end};
    
    %countyNormalizedDamages = countyDamages ./ countyDamages(:,1);
    
    countyInterpolatedDamages = 100*(countyDamages(:,1)' - interp1(countyDamageASLRs, countyDamages', [lowASLR, medianASLR, highASLR], 'linear')) ./ (countyDamages(:,1)');
     yNeg = countyInterpolatedDamages(2,:) - countyInterpolatedDamages(1,:);
    yPos = countyInterpolatedDamages(3,:) - countyInterpolatedDamages(2,:);
    
    
    for countyNum = 1:height(countyDamageTable)
        errorbar(countyInterpolatedDamages(2,countyNum),countyNum, yNeg(countyNum), yPos(countyNum), 'horizontal', 'o', 'LineWidth', 1.5, 'Color', colors(countyNum,:));
    
        %line([countyNum, countyNum], countyInterpolatedDamages([1 3],countyNum), 'LineWidth', 2);
    end

   
    
    %scatter(1:height(countyDamageTable), countyInterpolatedDamages(2,:))
    xlabel('Percentage of damage from anthropogenic SLR', 'fontsize', fontSize+4);
    
    %title('County Property Damage Differences');
    
    %set(gca,'ytick',[])
    set(gcf,'color','w');
    
    %set(gca, 'XDir','reverse')
    yticks(1:height(countyDamageTable))
   % yticklabels(strrep(countyDamageTable.PLACE_NAME(:), ' County', ''))
    ax = gca;
    for countyNum = 1:height(countyDamageTable)
        ax.YTickLabel{countyNum} = sprintf('\\color[rgb]{%f,%f,%f}%s', colors(countyNum,:), strrep(countyDamageTable.PLACE_NAME{countyNum}, ' County', ''));
    end
    
    a = get(gca,'yticklabel');
    set(gca, 'YtickLabel', a, 'fontsize', fontSize);
    
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print(['figures/county_damage_uncertainty'], '-dpdf', '-painters');
    
    
end

