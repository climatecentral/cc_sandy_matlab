% This code fits and sampling from the central budget
% estimate to use as the budget summary, which then is combined with the
% semi-empirical estimate to create the summary "ensemble" results.
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

% Code steps:
%   1. Budget summary
%       a. Load in/define the budget data
%       b. Fit a skew-normal distribution to the central estimate
%       c. Sample from the fit distribution N*3 times
%   2. Semi-empirical summary
%       a. Load in semi-empirical samples
%       b. Sample each budget with large N
%       c. Pool all samples and compute mean/CIs
%   3. Combine/Aggregate both summaries into an "Ensemble"
%       a. Pool all samples and compute mean/CIs
%       b. Print out the results into a Table
%

rng(1) % fix the random seed to ensure results are reproducable
N=1e8; % Define a large value for N (here we take 100 million)
q=[0.05,0.5,0.95]; % quantiles in our study

%% Step 1: Budget Summary Estimate

% 1a. Define the budget data
budget=struct('global', ...
                struct('low', [4.05214305, 6.965883269, 10.3771095], ...
                        'central', [6.441345278, 9.833290559, 13.74297191], ...
                        'high', [8.72945414, 12.62843735, 17.05264391]), ...
              'nyc', ...
                struct('low', [3.256081723, 6.127651696, 9.50292894], ...
                        'central', [5.241840801, 8.91736844, 13.07434067], ...
                        'high', [6.675407114, 11.63729272, 17.01820709]), ...
               'quantiles',q);

% convert to standard deviations, function below
Z90=1.65; % factor for converting from 90% CI (5th-95th) to standard deviation, assuming a normal distribution

global_central_s05=Xtostd(budget.global.central(1),budget.global.central(2),Z90);
global_central_s95=Xtostd(budget.global.central(3),budget.global.central(2),Z90);

nyc_central_s05=Xtostd(budget.nyc.central(1),budget.nyc.central(2),Z90);
nyc_central_s95=Xtostd(budget.nyc.central(3),budget.nyc.central(2),Z90);

% 1b. optimize a skew-normal distribution for the budget summary

% define the quantiles to fit
gmsl_budget_q=budget.global.central;
nyc_budget_q=budget.nyc.central;

% choose options for the optimization
options = optimset('Display','iter','MaxFunEvals',3e2,'TolFun',1e-05);
kurt0=3; % define the fixed kurtosis (for a standard normal distribution)

% find the initial guesses
gmsl0=[budget.global.central(2),mean([global_central_s05,global_central_s95]),0.];
nyc0=[budget.nyc.central(2),mean([nyc_central_s05,nyc_central_s95]),0.];

% define the noise magnitude to add to each moment
mu_noise=0.5;
std_noise=0.5;
skew_noise=0.5;

% define the points along which to evaluate the density
buffer=5;
gmsl_pts=round(budget.global.central(1),1)-buffer:0.01:round(budget.global.central(3),1)+buffer;
nyc_pts=round(budget.nyc.central(1),1)-buffer:0.01:round(budget.nyc.central(3),1)+buffer;

% define the arrays before the loop
nstart=50; % number of times to re-optimize with different initial conditions
ntest=1e4; % number of samples to fit with
rfit_gmsl=NaN(nstart,ntest); rfit_nyc=NaN(nstart,ntest);
qfit_gmsl=NaN(nstart,length(q)); qfit_nyc=NaN(nstart,length(q));
qerror_gmsl=NaN(nstart,1); qerror_nyc=NaN(nstart,1);
fdensity_gmsl=NaN(nstart,length(gmsl_pts)); fdensity_nyc=NaN(nstart,length(nyc_pts));
pfit_gmsl=NaN(nstart,3); pfit_nyc=NaN(nstart,3);

% loop to optimize for the best-fit LWS sample distribution
tic
for k=1:nstart
    
    % add noise to ensure we don't start in the same spot each iteration
    if k>1
        noise1 = unifrnd(-mu_noise,mu_noise);
        noise2 = unifrnd(0,std_noise);
        noise3 = unifrnd(-skew_noise,skew_noise);
        gmsl_in=[gmsl0(1)+noise1,gmsl0(2)+noise2,gmsl0(3)+noise3];
        nyc_in=[nyc0(1)+noise1,nyc0(2)+noise2,nyc0(3)+noise3];
    else
        gmsl_in=gmsl0;
        nyc_in=nyc0;
    end
    
    % fit the parameters
    fit_params = fminsearch(@fit_budget_gmsl,gmsl_in,options);
    pfit_gmsl(k,:)=squeeze(fit_params);
    clear fit_params
    
    fit_params = fminsearch(@fit_budget_nyc,nyc_in,options);
    pfit_nyc(k,:)=squeeze(fit_params);
    clear fit_params
    
    % sample a distribution from each optmized fit and find the quantiles
    rfit_gmsl(k,:)=squeeze(pearsrnd(pfit_gmsl(k,1),pfit_gmsl(k,2),pfit_gmsl(k,3),kurt0,ntest,1));
    qfit_gmsl(k,:)=squeeze(quantile(squeeze(rfit_gmsl(k,:)),q));
    
    rfit_nyc(k,:)=squeeze(pearsrnd(pfit_nyc(k,1),pfit_nyc(k,2),pfit_nyc(k,3),kurt0,ntest,1));
    qfit_nyc(k,:)=squeeze(quantile(squeeze(rfit_nyc(k,:)),q));
    
    % calculate the error with the original quantiles
    qerror_gmsl(k,1)=((gmsl_budget_q(1)-qfit_gmsl(k,1))^2 + (gmsl_budget_q(2)-qfit_gmsl(k,2))^2 + (gmsl_budget_q(3)-qfit_gmsl(k,3))^2);
    qerror_nyc(k,1)=((nyc_budget_q(1)-qfit_nyc(k,1))^2 + (nyc_budget_q(2)-qfit_nyc(k,2))^2 + (nyc_budget_q(3)-qfit_nyc(k,3))^2);
    
    % evaluate the density of the distribution
    fdensity_gmsl(k,:) = ksdensity(squeeze(rfit_gmsl(k,:)),gmsl_pts);
    fdensity_nyc(k,:) = ksdensity(squeeze(rfit_nyc(k,:)),nyc_pts);
    
    % print what iteration we are on
    k
    
end

% find the index of the best fit
[~,besti_gmsl]=min(qerror_gmsl)
[~,besti_nyc]=min(qerror_nyc)

% check the quantiles of the best fit
quantile(rfit_gmsl(besti_gmsl,:),q)
quantile(rfit_nyc(besti_nyc,:),q)

% time the optimization
toc

% plot the optimization
close all
h1=figure(1);
hold on
    plot(gmsl_pts,fdensity_gmsl(besti_gmsl,:),'k','LineWidth',2)
    plot(nyc_pts,fdensity_nyc(besti_nyc,:),'r','LineWidth',2)
    grid
    xlabel('GMSL/NYCSL in 2012 (cm)')
    ylabel('Density')
    title('Best-fit sample from GMSL/NYC Central Budgets (cm)')
    legend('GMSL','NYC')
    axis('tight')
hold off
set(h1,'Position', [500 150 700 500]);

% 1c. sample from the fit budget distributions
global_budgetsummary_samps = pearsrnd(pfit_gmsl(besti_gmsl,1),pfit_gmsl(besti_gmsl,2),pfit_gmsl(besti_gmsl,3),3,N*3,1);
nyc_budgetsummary_samps = pearsrnd(pfit_nyc(besti_nyc,1),pfit_nyc(besti_nyc,2),pfit_nyc(besti_nyc,3),3,N*3,1);

% double-check quantiles
gmsl_budget_centralsummary_q=quantile(global_budgetsummary_samps,q);
nyc_budget_centralsummary_q=quantile(nyc_budgetsummary_samps,q);

%% Step 2: Semi-empirical Summary Estimate

% 2a. Load in semi-empirical data
SEdatfile='./data/fig1/SEanalysis.mat';
SEdat=load(SEdatfile);
% define the time index to sample at
yrind=113; % yrind==113 in 2012

% 2b. Sample each budget with large number of samples, N

% Sample from each ASLR distribution
% HadCRUT4 Stable
stable_samps=datasample(SEdat.hadcrut.stable.summary.sldiff(:,yrind),N);
% HadCRUT4 Cooling
cooling_samps=datasample(SEdat.hadcrut.cooling.summary.sldiff(:,yrind),N);
% CMIP5 Counterfactual
cf_samps=datasample(SEdat.cmip5.cf.summary.sldiff(:,yrind),N);

% 2c. Aggregate and calculate mean and confidence intervals
aggregate_samps=vertcat(stable_samps,cooling_samps,cf_samps);
SE_mean=nanmean(aggregate_samps);
SE_50CI=quantile(aggregate_samps,q(2));
SE_05CI=quantile(aggregate_samps,q(1));
SE_95CI=quantile(aggregate_samps,q(3));

%% Step 3: Aggregate the Samples to Create the Ensemble

% 3a. Pool all samples and compute mean/CIs
% pool all samples together
global_summary_ALL_samps=vertcat(aggregate_samps,global_budgetsummary_samps);
% for NYC, scale the SE results with NYC scaling factor
nyc_scale_factor=0.91;
nyc_summary_ALL_samps=vertcat(aggregate_samps.*nyc_scale_factor,nyc_budgetsummary_samps);

% calculate mean and confidence intervals
global_ALL_mean=nanmean(global_summary_ALL_samps);
global_ALL_median=nanmedian(global_summary_ALL_samps);
global_ALL_05CI=quantile(global_summary_ALL_samps,q(1));
global_ALL_95CI=quantile(global_summary_ALL_samps,q(3));

nyc_ALL_mean=nanmean(nyc_summary_ALL_samps);
nyc_ALL_median=nanmedian(nyc_summary_ALL_samps);
nyc_ALL_05CI=quantile(nyc_summary_ALL_samps,q(1));
nyc_ALL_95CI=quantile(nyc_summary_ALL_samps,q(3));

% Look at Total semi-empirical Results pooled
cmipTOTAL_SAMPS=datasample(SEdat.cmip5.historical.summary.sl(:,yrind),N);
hadcrutTOTAL_SAMPS=datasample(SEdat.hadcrut.historical.summary.sl(:,yrind),N);
pool_se_gmsl=vertcat(cmipTOTAL_SAMPS,hadcrutTOTAL_SAMPS);
pool_se_gmsl_q=quantile(pool_se_gmsl,q);

%% 3b. Print to a Table

% generate the column names
summary_table_columns={'Scenario','GMSL50','GMSL5', 'GMSL95','NYC50','NYC5', 'NYC95'};

% generate the rows
summary_table_rows={'Central Budget (3/6)','Semi-Empirical (1/3+1/3+1/3)','Ensemble (1/6+1/6+1/6+3/6)'}';
Summary_table=table(summary_table_rows, ...
    [round(gmsl_budget_centralsummary_q(2),1); round(SE_50CI,1); round(global_ALL_median,1)], ...
    [round(gmsl_budget_centralsummary_q(1),1); round(SE_05CI,1); round(global_ALL_05CI,1)], ...
    [round(gmsl_budget_centralsummary_q(3),1); round(SE_95CI,1); round(global_ALL_95CI,1)], ...
    [round(nyc_budget_centralsummary_q(2),1); round(SE_50CI*nyc_scale_factor,1); round(nyc_ALL_median,1)], ...
    [round(nyc_budget_centralsummary_q(1),1); round(SE_05CI*nyc_scale_factor,1); round(nyc_ALL_05CI,1)], ...
    [round(nyc_budget_centralsummary_q(3),1); round(SE_95CI*nyc_scale_factor,1); round(nyc_ALL_95CI,1)], ...
    'VariableNames',summary_table_columns)
% writetable(Summary_table,'../Summary_table.xlsx')


%% Function Library

% sigma = (X - mean)/Z, and assuming
% median==mean for assumed normal distribution
function stdev = Xtostd(X,mean,Z)
    stdev = abs((X-mean)./Z);
end

% X = mean - z*sigma, and assuming
% median==mean for assumed normal distribution
function X = stdtoX_minus(sigma,mean,Z)
    X = mean - (Z.*sigma);
end

% X = mean + z*sigma, and assuming
% median==mean for assumed normal distribution
function X = stdtoX_plus(sigma,mean,Z)
    X = mean + (Z.*sigma);
end

% define the budget gmsl error function to optimize
function q_error=fit_budget_gmsl(m)
    % define the true quantile values we are seeking
    q_vals=[6.4, 9.8, 13.7];
    p05=q_vals(1);
    p50=q_vals(2);
    p95=q_vals(3);
    % set the constants for the search
    q=[0.05,0.5,0.95];
    n=1e4;
    kurt=3;

    % sample the distribution
    r=pearsrnd(m(1),m(2),m(3),kurt,n,1);
    
    % calculate the quantiles
    qstep=quantile(r,q);
    
    % calculate the difference
    q_error=(p05-qstep(1))^2 + (p50-qstep(2))^2 + (p95-qstep(3))^2;
    
end

% define the budget nyc error function to optimize
function q_error=fit_budget_nyc(m)
    % define the true quantile values we are seeking
    q_vals=[5.2, 8.9, 13.1];
    p05=q_vals(1);
    p50=q_vals(2);
    p95=q_vals(3);
    % set the constants for the search
    q=[0.05,0.5,0.95];
    n=1e5;
    kurt=3;

    % sample the distribution
    r=pearsrnd(m(1),m(2),m(3),kurt,n,1);
    
    % calculate the quantiles
    qstep=quantile(r,q);
    
    % calculate the difference
    q_error=(p05-qstep(1))^2 + (p50-qstep(2))^2 + (p95-qstep(3))^2;
    
end