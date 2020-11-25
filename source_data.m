% This code organizes and saves out part of the source data for publication.
%
% ---
% Code supporting Strauss et al. (2020) submitted to Nature
% Communications. If you use this code, please cite this study as:
%
%   B. H. Strauss, P. Orton, K. Bittermann, M. K. Buchanan, D. M. Gilford,
%   R. E. Kopp, S. Kulp, C. Massey, H. de Moel, S. Vinogradov, 2020: 
%   Economic Damages from Hurricane Sandy Attributable to Sea Level Rise 
%   Caused by Anthropogenic Climate Change. Nature Communications. (under
%   review, Nov. 2020)
%
% Input data and analyses are archived at [TK].
%
% Code credits:
%   -> Ben Strauss & Bob Kopp (project conception and development)
%   -> Klaus Bittermann & Bob Kopp (semi-empirical modeling)
%   -> Daniel Gilford (code development and maitenance)
% If you have any questions or comments, please contact Daniel Gilford at
% dgilford@climatecentral.org
% ---
% 

%% Setup
close all
clear

% define the save path
savepath='./data/SOURCE/';
q=[0.05,0.5,0.95];
nyf=0.91;
yrind=113;

%% Table 1

% Source file savename
T1_savename=fullfile(savepath,'Table 1/ASLR_estimates.xlsx');

% Define the Budget ASLR Data
budget=struct('global', ...
                struct('low', [4.05214305, 6.965883269, 10.3771095], ...
                        'central', [6.441345278, 9.833290559, 13.74297191], ...
                        'high', [8.72945414, 12.62843735, 17.05264391]), ...
              'nyc', ...
                struct('low', [3.256081723, 6.127651696, 9.50292894], ...
                        'central', [5.241840801, 8.91736844, 13.07434067], ...
                        'high', [6.675407114, 11.63729272, 17.01820709]));
                    
% load in the Semi-empirical data
SEdat=load('./data/fig1/SEanalysis.mat');

% load in the summary/ensemble data
summarydat=load('./data/source_data/summary_samps.mat');

% Columns
T1_col={'Scenario','Global50th','Global5th','Global95th','NY50th','NY5th','NY95th'};
% Rows
T1_row={'Budget Low Attribution','Budget Central Attribution','Budget High Attribution', ...
    'SE-Modeling Stable','SE-Modeling Cooling','SE-Modeling CMIP5', 'SE-Modeling Ensemble',...
    'Total Ensemble'};
% Build Table
T1=table(T1_row', ...
    [squeeze(budget.global.low(2));squeeze(budget.global.central(2));squeeze(budget.global.high(2));squeeze(SEdat.hadcrut.stable.summary.qdiff(2,yrind));squeeze(SEdat.hadcrut.cooling.summary.qdiff(2,yrind));squeeze(SEdat.cmip5.cf.summary.qdiff(2,yrind));squeeze(summarydat.SE_50CI);squeeze(summarydat.global_ALL_median)], ...
    [squeeze(budget.global.low(1));squeeze(budget.global.central(1));squeeze(budget.global.high(1));squeeze(SEdat.hadcrut.stable.summary.qdiff(1,yrind));squeeze(SEdat.hadcrut.cooling.summary.qdiff(1,yrind));squeeze(SEdat.cmip5.cf.summary.qdiff(1,yrind));squeeze(summarydat.SE_05CI);squeeze(summarydat.global_ALL_05CI)], ...
    [squeeze(budget.global.low(3));squeeze(budget.global.central(3));squeeze(budget.global.high(3));squeeze(SEdat.hadcrut.stable.summary.qdiff(3,yrind));squeeze(SEdat.hadcrut.cooling.summary.qdiff(3,yrind));squeeze(SEdat.cmip5.cf.summary.qdiff(3,yrind));squeeze(summarydat.SE_95CI);squeeze(summarydat.global_ALL_95CI)], ...
    [squeeze(budget.nyc.low(2));squeeze(budget.nyc.central(2));squeeze(budget.nyc.high(2));squeeze(SEdat.hadcrut.stable.summary.qdiff(2,yrind)*nyf);squeeze(SEdat.hadcrut.cooling.summary.qdiff(2,yrind)*nyf);squeeze(SEdat.cmip5.cf.summary.qdiff(2,yrind)*nyf);squeeze(summarydat.SE_50CI*nyf);squeeze(summarydat.nyc_ALL_median)], ...
    [squeeze(budget.nyc.low(1));squeeze(budget.nyc.central(1));squeeze(budget.nyc.high(1));squeeze(SEdat.hadcrut.stable.summary.qdiff(1,yrind)*nyf);squeeze(SEdat.hadcrut.cooling.summary.qdiff(1,yrind)*nyf);squeeze(SEdat.cmip5.cf.summary.qdiff(1,yrind)*nyf);squeeze(summarydat.SE_05CI*nyf);squeeze(summarydat.nyc_ALL_05CI)], ...
    [squeeze(budget.nyc.low(3));squeeze(budget.nyc.central(3));squeeze(budget.nyc.high(3));squeeze(SEdat.hadcrut.stable.summary.qdiff(3,yrind)*nyf);squeeze(SEdat.hadcrut.cooling.summary.qdiff(3,yrind)*nyf);squeeze(SEdat.cmip5.cf.summary.qdiff(3,yrind)*nyf);squeeze(summarydat.SE_95CI*nyf);squeeze(summarydat.nyc_ALL_95CI)], ...
    'VariableNames',T1_col);
writetable(T1,T1_savename);  


%% Figure 1

% load in the Figure data file
fig1dat=load('./data/source_data/fig1_data.mat');
yeargrid=fig1dat.fig1dat.axes.time;
nt=length(yeargrid)

% Source file savenames
TFig1a_savename=fullfile(savepath,'Fig 1/Fig1a_Timeseries.xlsx');
TFig1b_savename=fullfile(savepath,'Fig 1/Fig1b_Timeseries.xlsx');
TFig1c_savename=fullfile(savepath,'Fig 1/Fig1c_Timeseries.xlsx');

% ---------- Figure 1a ---------- %

% Fill Dangendorf Data
d19=ones(nt,1)*NaN;
d19_05=d19; d19_50=d19; d19_95=d19;
d19_05(end)=fig1dat.fig1dat.fig1a_historical.dangendorf2019(2);
d19_50(end)=fig1dat.fig1dat.fig1a_historical.dangendorf2019(1);
d19_95(end)=fig1dat.fig1dat.fig1a_historical.dangendorf2019(3);
% Fill Budget Data
bdg=ones(nt,1)*NaN;
bdg_05=d19; bdg_50=d19; bdg_95=d19;
bdg_05(end)=fig1dat.fig1dat.fig1a_historical.budget(1);
bdg_50(end)=fig1dat.fig1dat.fig1a_historical.budget(2);
bdg_95(end)=fig1dat.fig1dat.fig1a_historical.budget(3);

% Columns
TFig1a_col={'Time', ...
    'HadCRUT4_Hist_CI50','HadCRUT4_Hist_CI05','HadCRUT4_Hist_CI95', ...
    'CMIP5_Hist_CI50','CMIP5_Hist_CI05','CMIP5_Hist_CI95', ...
    'Dangendorf2019_CI50','Dangendorf2019_CI05','Dangendorf2019_CI95', ...
    'Budget_Hist_CI05','Budget_Hist_CI50','Budget_Hist_CI95'};

% Build Table
TFig1a=table(yeargrid', ...
    squeeze(fig1dat.fig1dat.fig1a_historical.hadcrut(2,:))',squeeze(fig1dat.fig1dat.fig1a_historical.hadcrut(1,:))',squeeze(fig1dat.fig1dat.fig1a_historical.hadcrut(3,:))', ...
    squeeze(fig1dat.fig1dat.fig1a_historical.cmip5(2,:))',squeeze(fig1dat.fig1dat.fig1a_historical.cmip5(1,:))',squeeze(fig1dat.fig1dat.fig1a_historical.cmip5(3,:))', ...
    squeeze(d19_50),squeeze(d19_05),squeeze(d19_95), ...
    squeeze(bdg_50),squeeze(bdg_05),squeeze(bdg_95), ...
    'VariableNames',TFig1a_col);
writetable(TFig1a,TFig1a_savename);

% ---------- Figure 1b ---------- %

% Columns
TFig1b_col={'Time', ...
    'HadCRUT4_Stable_CI50','HadCRUT4_Stable_CI05','HadCRUT4_Stable_CI95', ...
    'HadCRUT4_Cooling_CI50','HadCRUT4_Cooling_CI05','HadCRUT4_Cooling_CI95', ...
    'CMIP5_CF_CI50','CMIP5_CF_CI05','CMIP5_CF_CI95'};

% Build Table
TFig1b=table(yeargrid', ...
    squeeze(fig1dat.fig1dat.fig1b_counterfactual.hadcrut_stable(2,:))',squeeze(fig1dat.fig1dat.fig1b_counterfactual.hadcrut_stable(1,:))',squeeze(fig1dat.fig1dat.fig1b_counterfactual.hadcrut_stable(3,:))', ...
    squeeze(fig1dat.fig1dat.fig1b_counterfactual.hadcrut_cooling(2,:))',squeeze(fig1dat.fig1dat.fig1b_counterfactual.hadcrut_cooling(1,:))',squeeze(fig1dat.fig1dat.fig1b_counterfactual.hadcrut_cooling(3,:))', ...
    squeeze(fig1dat.fig1dat.fig1b_counterfactual.cmip5_cf(2,:))',squeeze(fig1dat.fig1dat.fig1b_counterfactual.cmip5_cf(1,:))',squeeze(fig1dat.fig1dat.fig1b_counterfactual.cmip5_cf(3,:))', ...
    'VariableNames',TFig1b_col);
writetable(TFig1b,TFig1b_savename);

% ---------- Figure 1c ---------- %

% Fill Budget Data
bdg_aslr=ones(nt,1)*NaN;
bdg_aslr_05=d19; bdg_aslr_50=d19; bdg_aslr_95=d19;
bdg_aslr_05(end)=fig1dat.fig1dat.fig1c_differences.budget(1);
bdg_aslr_50(end)=fig1dat.fig1dat.fig1c_differences.budget(2);
bdg_aslr_95(end)=fig1dat.fig1dat.fig1c_differences.budget(3);

% Columns
TFig1c_col={'Time', ...
    'HadCRUT4_Stable_ASLR_CI50','HadCRUT4_Stable_ASLR_CI05','HadCRUT4_Stable_ASLR_CI95', ...
    'HadCRUT4_Cooling_ASLR_CI50','HadCRUT4_Cooling_ASLR_CI05','HadCRUT4_Cooling_ASLR_CI95', ...
    'CMIP5_CF_ASLR_CI50','CMIP5_CF_ASLR_CI05','CMIP5_CF_ASLR_CI95', ...
    'Budget_Hist_ASLR_CI05','Budget_Hist_ASLR_CI50','Budget_Hist_ASLR_CI95'};

% Build Table
TFig1c=table(yeargrid', ...
    squeeze(fig1dat.fig1dat.fig1c_differences.hadcrut_stable_aslr(2,:))',squeeze(fig1dat.fig1dat.fig1c_differences.hadcrut_stable_aslr(1,:))',squeeze(fig1dat.fig1dat.fig1c_differences.hadcrut_stable_aslr(3,:))', ...
    squeeze(fig1dat.fig1dat.fig1c_differences.hadcrut_cooling_aslr(2,:))',squeeze(fig1dat.fig1dat.fig1c_differences.hadcrut_cooling_aslr(1,:))',squeeze(fig1dat.fig1dat.fig1c_differences.hadcrut_cooling_aslr(3,:))', ...
    squeeze(fig1dat.fig1dat.fig1c_differences.cmip5_cf_aslr(2,:))',squeeze(fig1dat.fig1dat.fig1c_differences.cmip5_cf_aslr(1,:))',squeeze(fig1dat.fig1dat.fig1c_differences.cmip5_cf_aslr(3,:))', ...
    squeeze(bdg_aslr_50),squeeze(bdg_aslr_05),squeeze(bdg_aslr_95), ...
    'VariableNames',TFig1c_col);
writetable(TFig1c,TFig1c_savename);


%% Table S1



%% Table S2

%% Table S4

%% Table S5

