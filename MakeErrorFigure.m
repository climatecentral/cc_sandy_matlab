fig = figure('Units','inches','position', [0, 0, 9, 5]);
set(fig,'PaperPosition', [0 0 9 9]);
set(fig,'Color','white');

base_ax = axes('Visible','off');
redblue = flipud(dlmread('redblue_256_rgb.txt')/255);
thermal = @(x) cmocean('thermal',x);

if ~exist('warpedError')
    warpedError = EvaluateSandyError( sandyErrorTable, sandyWarpedModelCompleteData, 3+6 );
end
figureRes=300;

lonlim = [-75.0 ,-71.8];
latlim = [38.8,42.0];


R = georasterref('Latlim',latlim,'Lonlim',lonlim,'RasterSize',[figureRes figureRes],'ColumnsStartFrom','south');
[figureLon, figureLat] = meshgrid(linspace(lonlim(1), lonlim(2), figureRes), linspace(latlim(1), latlim(2), figureRes));
states_hires = shaperead('E:/data/slr1/ss2/lidar/shoreline/shoreline_simple_noz.shp', 'UseGeoCoords', true);

states = shaperead('usastatehi', 'UseGeoCoords', true);

%if ~exist('figureRasterNonWarped')
    figureRasterNonWarped = modelFloodInterpolant{11}(figureLon,figureLat);
%end

%if ~exist('figureRasterWarped')
    %figureRasterWarped = modelFloodInterpolant{12}(figureLon,figureLat);
%end

%noDataZone = figureRasterWarped == figureRasterNonWarped;
%figureRasterNonWarped(noDataZone) = nan;
%figureRasterWarped(noDataZone) = nan;

%sp2 = subplot(1,3,1);
% 
% ax_2 = usamap(latlim,lonlim);
% set(ax_2, 'Visible', 'off')
% meshm(figureRasterNonWarped, R,[]);
% geoshow(ax_2, states_hires, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha',1, 'EdgeAlpha', 0)
% caxis([1 4])
% colormap(sp2, thermal(100));
% title({'Modeled', 'Max Water Elevation'});
% %geoshow(ax_1, states, 'FaceColor', [1 1 1], 'FaceAlpha',0)
% gridm off; mlabel off; plabel off
% %scatterm(sandyErrorTable.LAT_DEG_NORTH, sandyErrorTable.LON_DEG_EAST, 20, sandyErrorTable.OBS_METERS_NAVD88,'filled','o','MarkerEdgeColor','black','LineWidth',1);
% 
% 
% 
% sp1 = subplot(2,2,2);
% 
% 
% ax_1 = usamap(latlim,lonlim);
% set(ax_1, 'Visible', 'off')
% meshm(figureRasterNonWarped, R,[]);
% geoshow(ax_1, states_hires, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha',1, 'EdgeAlpha', 0)
% scatterm(sandyErrorTable.LAT_DEG_NORTH, sandyErrorTable.LON_DEG_EAST, 20, sandyErrorTable.OBS_METERS_NAVD88,'filled','o','MarkerEdgeColor','black','LineWidth',1);
% title({'Observed + Modeled', 'Max Water Elevation'});
% caxis([1 4])
% colormap(sp1, thermal(100));
% gridm off; mlabel off; plabel off




sp3 = subtightplot(1,2,1);



ax_3 = usamap(latlim,lonlim);
set(ax_3, 'Visible', 'off')

%meshm(figureRasterWarped, R,[]);
%caxis([1 5])
geoshow(ax_3, states, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha',1, 'EdgeColor', [0.8 0.8 0.8])
title({'Pre-Correction'});
gridm off; mlabel off; plabel off
scatterm(sandyErrorTable.LAT_DEG_NORTH, sandyErrorTable.LON_DEG_EAST, 20, sandyErrorTable.CC_MODEL_06_MINUS_OBS,'filled','o','MarkerEdgeColor','black','LineWidth',1);
caxis([-1.5 1.5])
colormap(sp3, redblue)

sp4 = subtightplot(1,2,2);
sp4.Position(1) = sp4.Position(1)-0.05;
ax_4 = usamap(latlim,lonlim);
set(ax_4, 'Visible', 'off')
title({'Post-Correction'});
geoshow(ax_4, states, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha',1, 'EdgeColor', [0.8 0.8 0.8])
gridm off; mlabel off; plabel off
scatterm(sandyErrorTable.LAT_DEG_NORTH, sandyErrorTable.LON_DEG_EAST, 20, warpedError(:,3),'filled','o','MarkerEdgeColor','black','LineWidth',1);
caxis([-1.5 1.5])

colormap(sp4, redblue)


hp4 = sp4.Position;%get(sp4,'Position');
cb = colorbar('Position', [hp4(1)+hp4(3)-0.01  hp4(2)  0.05  hp4(2)+hp4(3)*2]);
cb.TickLabels = {'-150', '-100', '-50', '0', '50', '100', '150'};

export_fig('images/floodmap_legend','-opengl','-RGB','-r300','-nocrop',gcf);

    

    
      set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print(['figures/error_fig'], '-dpdf', '-painters');
    sandyErrorTable = removevars(sandyErrorTable, {'MODEL_01_MAXELEV_NAVD88','CC_MODEL_01_MINUS_OBS','MODEL_02_MAXELEV_NAVD88','CC_MODEL_02_MINUS_OBS','MODEL_03_MAXELEV_NAVD88','CC_MODEL_03_MINUS_OBS','MODEL_04_MAXELEV_NAVD88','CC_MODEL_04_MINUS_OBS','MODEL_05_MAXELEV_NAVD88','CC_MODEL_05_MINUS_OBS','MODEL_06_MAXELEV_NAVD88','CC_MODEL_06_MINUS_OBS','MODEL_07_MAXELEV_NAVD88','CC_MODEL_07_MINUS_OBS','MODEL_08_MAXELEV_NAVD88','CC_MODEL_08_MINUS_OBS','MODEL_09_MAXELEV_NAVD88','CC_MODEL_09_MINUS_OBS'});
sandyErrorTable = removevars(sandyErrorTable, 'Var15');

    sandyErrorTable.PreCorrectionError = nowarpError;
    sandyErrorTable.PostCorrectionError = warpedError;
 