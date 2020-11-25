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
nyf_noLWS=0.883761834;
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
TFig1a_savename=fullfile(savepath,'MS Fig 1/Fig1a_Timeseries.xlsx');
TFig1b_savename=fullfile(savepath,'MS Fig 1/Fig1b_Timeseries.xlsx');
TFig1c_savename=fullfile(savepath,'MS Fig 1/Fig1c_Timeseries.xlsx');

% ---------- Figure 1a ---------- %

% Fill Dangendorf Data
d19=ones(nt,1)*NaN;
d19_05=d19; d19_50=d19; d19_95=d19;
d19_05(end)=fig1dat.fig1dat.fig1a_historical.dangendorf2019(2);
d19_50(end)=fig1dat.fig1dat.fig1a_historical.dangendorf2019(1);
d19_95(end)=fig1dat.fig1dat.fig1a_historical.dangendorf2019(3);
% Fill Budget Data
bdg=ones(nt,1)*NaN;
bdg_05=bdg; bdg_50=bdg; bdg_95=bdg;
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

% Source file savename
TS1_savename=fullfile(savepath,'SI Table 1/Total_GMSL.xlsx');

% Define the LWS terms (from Budget)
lws_global=[0.751699404,-0.110000000,0.531699404]; lws_nyc=[0.574222152,1.092179865,1.610137578];
% Calculate the SE modeling + LWS
% HadCRUT4
hadcrut_noLWS=SEdat.hadcrut.historical.summary.q(:,yrind);
hadcrut_LWS(2)=hadcrut_noLWS(2)+lws_global(2);
hadcrut_LWS(1)=hadcrut_LWS(2)-sqrt((lws_global(2)-lws_global(1))^2+(hadcrut_noLWS(2)-hadcrut_noLWS(1))^2);
hadcrut_LWS(3)=hadcrut_LWS(2)+sqrt((lws_global(2)-lws_global(3))^2+(hadcrut_noLWS(2)-hadcrut_noLWS(3))^2);
% CMIP5
cmip5_noLWS=SEdat.cmip5.historical.summary.q(:,yrind);
cmip5_LWS(2)=cmip5_noLWS(2)+lws_global(2);
cmip5_LWS(1)=cmip5_LWS(2)-sqrt((lws_global(2)-lws_global(1))^2+(cmip5_noLWS(2)-cmip5_noLWS(1))^2);
cmip5_LWS(3)=cmip5_LWS(2)+sqrt((lws_global(2)-lws_global(3))^2+(cmip5_noLWS(2)-cmip5_noLWS(3))^2);

% Columns
TS1_col={'Scenario','Global50th','Global5th','Global95th'};
% Rows
TS1_row={'Observed','Budget', ...
    'SE-Modeling HadCRUT4+LWS','SE-Modeling CMIP5+LWS'};

% Build Table
TS1=table(TS1_row', ...
    [squeeze(fig1dat.fig1dat.fig1a_historical.dangendorf2019(1));squeeze(fig1dat.fig1dat.fig1a_historical.budget(2));hadcrut_LWS(2);cmip5_LWS(2)], ...
    [squeeze(fig1dat.fig1dat.fig1a_historical.dangendorf2019(2));squeeze(fig1dat.fig1dat.fig1a_historical.budget(1));hadcrut_LWS(1);cmip5_LWS(1)], ...
    [squeeze(fig1dat.fig1dat.fig1a_historical.dangendorf2019(3));squeeze(fig1dat.fig1dat.fig1a_historical.budget(3));hadcrut_LWS(3);cmip5_LWS(3)], ...
    'VariableNames',TS1_col);
writetable(TS1,TS1_savename);

%% Table S2

% Source file savename
TS2_savename=fullfile(savepath,'SI Table 2/Total_NYC.xlsx');

% Define the NOAA NY observations (copied from observed.m)
ny_observed=[30.4941,33.7488,37.0035];
% Define the GIA term (from Kopp 2013)
gia=[10.88490,14.69000,18.49510];
% Calculate Observed excluding GIA
ny_noGIA(2)=ny_observed(2)-gia(2);
ny_noGIA(1)=ny_noGIA(2)-sqrt((gia(2)-gia(1))^2+(ny_observed(2)-ny_observed(1))^2);
ny_noGIA(3)=ny_noGIA(2)+sqrt((gia(3)-gia(2))^2+(ny_observed(2)-ny_observed(3))^2);

% Define the NY total Budget Data
ny_total_bdg=[10.95605287,16.67771649,22.81673570];

% Calculate the SE modeling + LWS
% HadCRUT4
NYhadcrut_noLWS=SEdat.hadcrut.historical.summary.q(:,yrind)*nyf_noLWS;
NYhadcrut_LWS(2)=NYhadcrut_noLWS(2)+lws_nyc(2);
NYhadcrut_LWS(1)=NYhadcrut_LWS(2)-sqrt((lws_nyc(2)-lws_nyc(1))^2+(NYhadcrut_noLWS(2)-NYhadcrut_noLWS(1))^2);
NYhadcrut_LWS(3)=NYhadcrut_LWS(2)+sqrt((lws_nyc(2)-lws_nyc(3))^2+(NYhadcrut_noLWS(2)-NYhadcrut_noLWS(3))^2);
% CMIP5
NYcmip5_noLWS=SEdat.cmip5.historical.summary.q(:,yrind)*nyf_noLWS;
NYcmip5_LWS(2)=NYcmip5_noLWS(2)+lws_nyc(2);
NYcmip5_LWS(1)=NYcmip5_LWS(2)-sqrt((lws_nyc(2)-lws_nyc(1))^2+(NYcmip5_noLWS(2)-NYcmip5_noLWS(1))^2);
NYcmip5_LWS(3)=NYcmip5_LWS(2)+sqrt((lws_nyc(2)-lws_nyc(3))^2+(NYcmip5_noLWS(2)-NYcmip5_noLWS(3))^2);

% Columns
TS2_col={'Scenario','NY50th','NY5th','NY95th'};
% Rows
TS2_row={'Observed_noGIA','Budget', ...
    'SE-Modeling HadCRUT4+LWS','SE-Modeling CMIP5+LWS'};

% Build Table
TS2=table(TS2_row', ...
    [squeeze(ny_noGIA(2));squeeze(ny_total_bdg(2));squeeze(NYhadcrut_LWS(2));squeeze(NYcmip5_LWS(2))], ...
    [squeeze(ny_noGIA(1));squeeze(ny_total_bdg(1));squeeze(NYhadcrut_LWS(1));squeeze(NYcmip5_LWS(1))], ...
    [squeeze(ny_noGIA(3));squeeze(ny_total_bdg(3));squeeze(NYhadcrut_LWS(3));squeeze(NYcmip5_LWS(3))], ...
    'VariableNames',TS2_col)
writetable(TS2,TS2_savename);

%% Table S3

% This is copied directly from the budget table and added to a speadsheet;
% see Supplementary Table S3.

%% Table S4

% Source file savename
TS4_savename=fullfile(savepath,'SI Table 4/Global_SEmodel_Estimates.xlsx');

% Columns
TS4_col={'Scenario','Summary50th','Summary5th','Summary95th','Mann50th','Mann5th','Mann95th','Marcott50th','Marcott5th','Marcott95th'};
% Rows
TS4_row={'Observed',...
    'CMIP5 Hist.','CMIP5 CF',...
    'HadCRUT4 Hist.','HadCRUT4 Stable','HadCRUT4 Cooling',...
    'CMIP5 ASLR','HadCRUT4 Stable ASLR','HadCRUT4 Cooling ASLR'};

% Build Table
TS4=table(TS4_row', ...
    [squeeze(fig1dat.fig1dat.fig1a_historical.dangendorf2019(1));squeeze(SEdat.cmip5.historical.summary.q(2,yrind)); squeeze(SEdat.cmip5.cf.summary.q(2,yrind));squeeze(SEdat.hadcrut.historical.summary.q(2,yrind)); squeeze(SEdat.hadcrut.stable.summary.q(2,yrind)); squeeze(SEdat.hadcrut.cooling.summary.q(2,yrind));squeeze(SEdat.cmip5.cf.summary.qdiff(2,yrind)); squeeze(SEdat.hadcrut.stable.summary.qdiff(2,yrind)); squeeze(SEdat.hadcrut.cooling.summary.qdiff(2,yrind))], ...
    [squeeze(fig1dat.fig1dat.fig1a_historical.dangendorf2019(2));squeeze(SEdat.cmip5.historical.summary.q(1,yrind)); squeeze(SEdat.cmip5.cf.summary.q(1,yrind));squeeze(SEdat.hadcrut.historical.summary.q(1,yrind)); squeeze(SEdat.hadcrut.stable.summary.q(1,yrind)); squeeze(SEdat.hadcrut.cooling.summary.q(1,yrind));squeeze(SEdat.cmip5.cf.summary.qdiff(1,yrind)); squeeze(SEdat.hadcrut.stable.summary.qdiff(1,yrind)); squeeze(SEdat.hadcrut.cooling.summary.qdiff(1,yrind))], ...
    [squeeze(fig1dat.fig1dat.fig1a_historical.dangendorf2019(3));squeeze(SEdat.cmip5.historical.summary.q(3,yrind)); squeeze(SEdat.cmip5.cf.summary.q(3,yrind));squeeze(SEdat.hadcrut.historical.summary.q(3,yrind)); squeeze(SEdat.hadcrut.stable.summary.q(3,yrind)); squeeze(SEdat.hadcrut.cooling.summary.q(3,yrind));squeeze(SEdat.cmip5.cf.summary.qdiff(3,yrind)); squeeze(SEdat.hadcrut.stable.summary.qdiff(3,yrind)); squeeze(SEdat.hadcrut.cooling.summary.qdiff(3,yrind))], ...
    [NaN;squeeze(SEdat.cmip5.historical.mann.q(2,yrind)); squeeze(SEdat.cmip5.cf.mann.q(2,yrind));squeeze(SEdat.hadcrut.historical.mann.q(2,yrind)); squeeze(SEdat.hadcrut.stable.mann.q(2,yrind)); squeeze(SEdat.hadcrut.cooling.mann.q(2,yrind));squeeze(SEdat.cmip5.cf.mann.qdiff(2,yrind)); squeeze(SEdat.hadcrut.stable.mann.qdiff(2,yrind)); squeeze(SEdat.hadcrut.cooling.mann.qdiff(2,yrind))], ...
    [NaN;squeeze(SEdat.cmip5.historical.mann.q(1,yrind)); squeeze(SEdat.cmip5.cf.mann.q(1,yrind));squeeze(SEdat.hadcrut.historical.mann.q(1,yrind)); squeeze(SEdat.hadcrut.stable.mann.q(1,yrind)); squeeze(SEdat.hadcrut.cooling.mann.q(1,yrind));squeeze(SEdat.cmip5.cf.mann.qdiff(1,yrind)); squeeze(SEdat.hadcrut.stable.mann.qdiff(1,yrind)); squeeze(SEdat.hadcrut.cooling.mann.qdiff(1,yrind))], ...
    [NaN;squeeze(SEdat.cmip5.historical.mann.q(3,yrind)); squeeze(SEdat.cmip5.cf.mann.q(3,yrind));squeeze(SEdat.hadcrut.historical.mann.q(3,yrind)); squeeze(SEdat.hadcrut.stable.mann.q(3,yrind)); squeeze(SEdat.hadcrut.cooling.mann.q(3,yrind));squeeze(SEdat.cmip5.cf.mann.qdiff(3,yrind)); squeeze(SEdat.hadcrut.stable.mann.qdiff(3,yrind)); squeeze(SEdat.hadcrut.cooling.mann.qdiff(3,yrind))], ...
    [NaN;squeeze(SEdat.cmip5.historical.marcott.q(2,yrind)); squeeze(SEdat.cmip5.cf.marcott.q(2,yrind));squeeze(SEdat.hadcrut.historical.marcott.q(2,yrind)); squeeze(SEdat.hadcrut.stable.marcott.q(2,yrind)); squeeze(SEdat.hadcrut.cooling.marcott.q(2,yrind));squeeze(SEdat.cmip5.cf.marcott.qdiff(2,yrind)); squeeze(SEdat.hadcrut.stable.marcott.qdiff(2,yrind)); squeeze(SEdat.hadcrut.cooling.marcott.qdiff(2,yrind))], ...
    [NaN;squeeze(SEdat.cmip5.historical.marcott.q(1,yrind)); squeeze(SEdat.cmip5.cf.marcott.q(1,yrind));squeeze(SEdat.hadcrut.historical.marcott.q(1,yrind)); squeeze(SEdat.hadcrut.stable.marcott.q(1,yrind)); squeeze(SEdat.hadcrut.cooling.marcott.q(1,yrind));squeeze(SEdat.cmip5.cf.marcott.qdiff(1,yrind)); squeeze(SEdat.hadcrut.stable.marcott.qdiff(1,yrind)); squeeze(SEdat.hadcrut.cooling.marcott.qdiff(1,yrind))], ...
    [NaN;squeeze(SEdat.cmip5.historical.marcott.q(3,yrind)); squeeze(SEdat.cmip5.cf.marcott.q(3,yrind));squeeze(SEdat.hadcrut.historical.marcott.q(3,yrind)); squeeze(SEdat.hadcrut.stable.marcott.q(3,yrind)); squeeze(SEdat.hadcrut.cooling.marcott.q(3,yrind));squeeze(SEdat.cmip5.cf.marcott.qdiff(3,yrind)); squeeze(SEdat.hadcrut.stable.marcott.qdiff(3,yrind)); squeeze(SEdat.hadcrut.cooling.marcott.qdiff(3,yrind))], ...
    'VariableNames',TS4_col)
writetable(TS4,TS4_savename);


%% Table S5

% Source file savename
TS5_savename=fullfile(savepath,'SI Table 5/CMIP5_IndividualModel_Estimates.xlsx');

% Reorganize the model names, drop out GFDL
dropname='GFDL-CM3';
nmod_noDROP=length(SEdat.cmip5.metadata);
model_names={''};
count=1;
for m=1:nmod_noDROP
    modname=SEdat.cmip5.metadata(m).model_names;
    if strcmp(modname,dropname)==1
        continue
    else
        model_names{count}=strcat(SEdat.cmip5.metadata(m).model_names,{' '},SEdat.cmip5.metadata(m).replicate_names);
        count=count+1;
    end
end
model_names;

% find the quantiles for each individual model
% SUMMARY
Hist_Summary_q=squeeze(quantile(SEdat.cmip5.historical.summary.slmodels(:,yrind,:),q,1));
CF_Summary_q=squeeze(quantile(SEdat.cmip5.cf.summary.slmodels(:,yrind,:),q,1));
ASLR_Summary_q=squeeze(quantile(SEdat.cmip5.cf.summary.sldiffmodels(:,yrind,:),q,1));
% MARCOTT
Hist_Marcott_q=squeeze(quantile(SEdat.cmip5.historical.marcott.slmodels(:,yrind,:),q,1));
CF_Marcott_q=squeeze(quantile(SEdat.cmip5.cf.marcott.slmodels(:,yrind,:),q,1));
ASLR_Marcott_q=squeeze(quantile(SEdat.cmip5.cf.marcott.sldiffmodels(:,yrind,:),q,1));
% MANN
Hist_Mann_q=squeeze(quantile(SEdat.cmip5.historical.mann.slmodels(:,yrind,:),q,1));
CF_Mann_q=squeeze(quantile(SEdat.cmip5.cf.mann.slmodels(:,yrind,:),q,1));
ASLR_Mann_q=squeeze(quantile(SEdat.cmip5.cf.mann.sldiffmodels(:,yrind,:),q,1));

% Columns
TS5_col={'Model_Ensemble', ...
    'Summary_Hist_50th','Summary_Hist_5th','Summary_Hist_95th', ...
    'Summary_CF_50th','Summary_CF_5th','Summary_CF_95th', ...
    'Summary_ASLR_50th','Summary_ASLR_5th','Summary_ASLR_95th', ...
    'Mann_Hist_50th','Mann_Hist_5th','Mann_Hist_95th', ...
    'Mann_CF_50th','Mann_CF_5th','Mann_CF_95th', ...
    'Mann_ASLR_50th','Mann_ASLR_5th','Mann_ASLR_95th', ...
    'Marcott_Hist_50th','Marcott_Hist_5th','Marcott_Hist_95th', ...
    'Marcott_CF_50th','Marcott_CF_5th','Marcott_CF_95th', ...
    'Marcott_ASLR_50th','Marcott_ASLR_5th','Marcott_ASLR_95th'};

% Rows
TS5_row=model_names;

% Build Table
TS5=table([TS5_row{:}]', ...
    squeeze(Hist_Summary_q(2,:)'),squeeze(Hist_Summary_q(1,:)'),squeeze(Hist_Summary_q(3,:)'), squeeze(CF_Summary_q(2,:)'),squeeze(CF_Summary_q(1,:)'),squeeze(CF_Summary_q(3,:)'), squeeze(ASLR_Summary_q(2,:)'),squeeze(ASLR_Summary_q(1,:)'),squeeze(ASLR_Summary_q(3,:)'), ...
    squeeze(Hist_Mann_q(2,:)'),squeeze(Hist_Mann_q(1,:)'),squeeze(Hist_Mann_q(3,:)'), squeeze(CF_Mann_q(2,:)'),squeeze(CF_Mann_q(1,:)'),squeeze(CF_Mann_q(3,:)'), squeeze(ASLR_Mann_q(2,:)'),squeeze(ASLR_Mann_q(1,:)'),squeeze(ASLR_Mann_q(3,:)'), ...
    squeeze(Hist_Marcott_q(2,:)'),squeeze(Hist_Marcott_q(1,:)'),squeeze(Hist_Marcott_q(3,:)'), squeeze(CF_Marcott_q(2,:)'),squeeze(CF_Marcott_q(1,:)'),squeeze(CF_Marcott_q(3,:)'), squeeze(ASLR_Marcott_q(2,:)'),squeeze(ASLR_Marcott_q(1,:)'),squeeze(ASLR_Marcott_q(3,:)'), ...
    'VariableNames',TS5_col);
writetable(TS5,TS5_savename);
