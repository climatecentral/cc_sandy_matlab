% This code plots Figure 1, showing the observed and semi-empirical
% SLR timeseries, alongside the global SLR budget percentiles.
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

% define the quantiles and timegrid
timegrid=1900:1:2012;
nt=length(timegrid);
q=[0.05, 0.50, 0.95];

%% Observed Data (Dangendorf et al. 2019)

% load Dangendorf et al. (2019) timeseries data
% Taken from:
%   Dangendorf, S., Hay, C., Calafat, F.M. et al. 
%   Persistent acceleration in global sea-level rise since the 1960s. 
%   Nat. Clim. Chang. 9, 705?710 (2019). https://doi.org/10.1038/s41558-019-0531-8
d19_dat=importdata('./data/fig1/Dangendorf2019_GMSL.txt');
d19_time=d19_dat.data(:,1);

% convert to cm and reframe around 0 in 1900
d19_mean_unshifted=d19_dat.data(:,2)./10;
d19_mean=d19_mean_unshifted-d19_mean_unshifted(1);
d19_sigma=d19_dat.data(:,3)./10;

% convert the data from mon to yr,mon format
nyear=floor(d19_time(end))-floor(d19_time(1))+1;
d19_timegrid=floor(d19_time(1)):1:floor(d19_time(end));
nmon=12;
d19_mean_ym=zeros(nyear,nmon);
d19_sigma_ym=zeros(nyear,nmon);
for y=1:nyear
    for m=1:nmon
       if (y==1) & (m==1)
            d19_mean_ym(y,m)=NaN;
        d19_sigma_ym(y,m)=NaN;
       else
        tind=(y-1)*12+m;
        d19_mean_ym(y,m)=d19_mean(tind-1);
        d19_sigma_ym(y,m)=d19_sigma(tind-1);
       end
    end
end

% get the annual average
d19_annual_mean=squeeze(nanmean(d19_mean_ym,2));
d19_annual_sigma=squeeze(nanmean(d19_sigma_ym,2));

% get the 90th confidence interval
Z90=1.65; Z95=1.96;
d19_90_upper = stdtoX_plus(d19_annual_sigma,d19_annual_mean,Z90);
d19_90_lower = stdtoX_minus(d19_annual_sigma,d19_annual_mean,Z90);

% plot for example by uncommenting
figure(999)
hold on
    %plot(hay_dat(:,1),hay_dat(:,2),'c','LineWidth',2)
    plot(d19_timegrid,d19_annual_mean,'b','LineWidth',2)
    plot(d19_timegrid,d19_90_upper,'b--')
    plot(d19_timegrid,d19_90_lower,'b--')
    legend('Annual Mean','95th CI','Location','Best')
    ylabel('Global SLR (cm)')
    xlabel('Year')
hold off

% 2012 global SLR from Dangendorf Trends (copied from observed.m)
d19_range=[17.9200   10.5280   25.3120]

%% Semi-empirical results

% Load semi-empirical model results analyzed in SEanalysis.m
sedatfile='./data/fig1/SEanalysis.mat';
SEdat=load(sedatfile);

% find the quantile timeseries
% HISTORICAL (1a)
hadcrut_hist=quantile(SEdat.hadcrut.historical.summary.sl(:,1:nt),q);
cmip_hist=quantile(SEdat.cmip5.historical.summary.sl(:,1:nt),q);
% COUNTERFACTUAL (1b)
hadcrut_stable=quantile(SEdat.hadcrut.stable.summary.sl(:,1:nt),q);
hadcrut_cooling=quantile(SEdat.hadcrut.cooling.summary.sl(:,1:nt),q);
cmip_counterfactual=quantile(SEdat.cmip5.cf.summary.sl(:,1:nt),q);
% DIFFERENCES/ASLR (1c)
hadcrut_stable_aslr=quantile(SEdat.hadcrut.stable.summary.sldiff(:,1:nt),q);
hadcrut_cooling_aslr=quantile(SEdat.hadcrut.cooling.summary.sldiff(:,1:nt),q);
cmip_counterfactual_aslr=quantile(SEdat.cmip5.cf.summary.sldiff(:,1:nt),q);


%% SLR Budget

% define Budget values to add to figures 1a+c (from budget analyses,
% copied from Supplementary Table 3)
budget=struct('central_gmsl', [12.4, 17.5, 23.1], ...
               'central_aslr', [6.441345278, 9.833290559, 13.74297191]);

%% Plot all data together (Figure 1)

% CAPTION: Fig. 1. Curves of semi-empirically modeled (a) historical and 
% (b) counterfactual GMSL changes from 1900?2012, and(c) their differences. 
% The differences estimate anthropogenic sea level rise; shaded regions 
% indicate 5th?95th percentileranges. Differences are taken between matching 
% scenarios: historical vs. natural CMIP5 simulations; and historical HadCRUT4-based 
% scenario vs. stable or cooling counterfactual temperature-based scenarios from 
% Kopp et al (2016).Curves  reflect  results  pooled  across  calibrations  for  
% both  temperature  reconstructions  (Mann  et  al,  2009;  Marcottet  al,  2013).  
% Also  included  are  the  medians  (points)  and  5th?95th  percentiles  (vertical 
% bars)  of  historical  (panel a)  and  central-case  attributable  (panel  c)  sea  
% level  rise  budgets  (Supplementary  Table  3).  Observed  GMSL  from Dangendorf 
% et al (2019) over 1900?2012 (median and 5th?95th percentiles) is also shown for 
% comparison (panel a; Supplementary Table 1).

% figure details
set(0,'defaultAxesFontSize',13)
faceAlpha=0.2;
figsizex=250;
figsizey=500;
set(gca, 'FontName', 'Arial')

% Figure 1a: Observed/Historical SLR
% Include HadCRUT4, CMIP5, and Observed (Dangendorf et al. 2019)
h1=figure(1);
set(gcf,'Color','w')
hold on
    X = [timegrid fliplr(timegrid)];
    
    % CMIP
    Y = [cmip_hist(1,:) fliplr(cmip_hist(3,:))];
    fill(X,Y,'m','EdgeAlpha',0,'FaceAlpha',faceAlpha,'HandleVisibility','off')
    plot(timegrid,cmip_hist(2,:),'color','m','LineWidth',3)
    
    % HadCRUT
    Y = [hadcrut_hist(1,:) fliplr(hadcrut_hist(3,:))];
    fill(X,Y,[1 0.5 0],'EdgeAlpha',0,'FaceAlpha',faceAlpha,'HandleVisibility','off')
    plot(timegrid,hadcrut_hist(2,:),'color',[1 0.5 0],'LineWidth',3)
    
    % Dangendorf et al. (2019)
    plot([2015 2015],[d19_range(2) d19_range(3)],'Color',[0.5 0.5 0.5],'LineWidth',5,'HandleVisibility','off')
    plot(2015,d19_range(1),'o','MarkerSize',12,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w')
    
    % Budget
    plot([2020 2020],[budget.central_gmsl(1) budget.central_gmsl(3)],'k','LineWidth',5,'HandleVisibility','off')
    plot(2020,budget.central_gmsl(2),'o','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','w')
    
    legend('CMIP5','HadCRUT4','Observed','Budget','Location','Northwest')
    
    xlim([1900,2020])
    ylim([-3,26])
    xticks([1900:20:2010])
    xlabel('Year')
    ylabel('Global Mean Sea Level Rise (cm)')
    title('a)   Historical Global Mean Sea Level Rise')

hold off
set(gcf,'Position',[100 100 figsizex figsizey])
% saveas(h1,'fig1a.pdf')

% Figure 1b: Counterfactual
% Include HadCRUT4 Cooling, Stable, CMIP5 CF
h2=figure(2);
set(gcf,'Color','w')
hold on
    X = [timegrid fliplr(timegrid)];
    
    % CMIP5 CF
    Y = [cmip_counterfactual(1,:) fliplr(cmip_counterfactual(3,:))];
    fill(X,Y,'m','EdgeAlpha',0,'FaceAlpha',faceAlpha,'HandleVisibility','off')
    plot(timegrid,cmip_counterfactual(2,:),'color','m','LineWidth',3,'HandleVisibility','off')
    
    % HadCRUT Stable
    Y = [hadcrut_stable(1,:) fliplr(hadcrut_stable(3,:))];
    fill(X,Y,'b','EdgeAlpha',0,'FaceAlpha',faceAlpha,'HandleVisibility','off')
    plot(timegrid,hadcrut_stable(2,:),'color','b','LineWidth',3)
    
    % HadCRUT Cooling
    Y = [hadcrut_cooling(1,:) fliplr(hadcrut_cooling(3,:))];
    fill(X,Y,[0,0.5,0],'EdgeAlpha',0,'FaceAlpha',faceAlpha,'HandleVisibility','off')
    plot(timegrid,hadcrut_cooling(2,:),'color',[0,0.5,0],'LineWidth',3)
    
    legend('Stable','Cooling','Location','Northwest')
    
    xlim([1900,2012])
    ylim([-3,26])
    xticks([1900:20:2020])
    xlabel('Year')
    ylabel('Global Mean Sea Level Rise (cm)')
    title('b)   Counterfactual')

hold off
set(gcf,'Position',[200 150 figsizex figsizey])
% saveas(h2,'fig1b.pdf')


% Figure 1c: Differences
% Include HadCRUT4 Cooling, Stable, CMIP5 CF, Budget
h3=figure(3);
set(gcf,'Color','w')
hold on
    X = [timegrid fliplr(timegrid)];
    
    % CMIP5 CF
    Y = [cmip_counterfactual_aslr(1,:) fliplr(cmip_counterfactual_aslr(3,:))];
    fill(X,Y,'m','EdgeAlpha',0,'FaceAlpha',faceAlpha,'HandleVisibility','off')
    plot(timegrid,cmip_counterfactual_aslr(2,:),'color','m','LineWidth',3)
    
    % HadCRUT Stable
    Y = [hadcrut_stable_aslr(1,:) fliplr(hadcrut_stable_aslr(3,:))];
    fill(X,Y,'b','EdgeAlpha',0,'FaceAlpha',faceAlpha,'HandleVisibility','off')
    plot(timegrid,hadcrut_stable_aslr(2,:),'color','b','LineWidth',3)
    
    % HadCRUT Cooling
    Y = [hadcrut_cooling_aslr(1,:) fliplr(hadcrut_cooling_aslr(3,:))];
    fill(X,Y,[0,0.5,0],'EdgeAlpha',0,'FaceAlpha',faceAlpha,'HandleVisibility','off')
    plot(timegrid,hadcrut_cooling_aslr(2,:),'color',[0,0.5,0],'LineWidth',3)
    
    % Budget
    plot([2015 2015],[budget.central_aslr(1) budget.central_aslr(3)],'k','LineWidth',5,'HandleVisibility','off')
    plot(2015,budget.central_aslr(2),'ko','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','w')
    
    xlim([1900,2015])
    ylim([-3,26])
    xticks([1900:20:2020])
    xlabel('Year')
    ylabel('Attributable Sea Level Rise (cm)')
    title('c)   Difference (ASLR)')

hold off
set(gcf,'Position',[300 200 figsizex figsizey])
% saveas(h3,'fig1c.pdf')

%% Save the data out for use in the source file

% save matlab structure
fig1dat=struct('axes', ...
                struct('time',timegrid, ...
                        'quantiles',q), ...
            'fig1a_historical', ...
                struct('hadcrut',hadcrut_hist, ...
                        'cmip5',cmip_hist, ...
                        'budget',budget.central_gmsl, ...
                        'dangendorf2019',d19_range), ...
            'fig1b_counterfactual', ...
                struct('hadcrut_cooling',hadcrut_cooling, ...
                        'hadcrut_stable',hadcrut_stable, ...
                        'cmip5_cf',cmip_counterfactual), ...
             'fig1c_differences', ...
                struct('hadcrut_cooling_aslr',hadcrut_cooling_aslr, ...
                        'hadcrut_stable_aslr',hadcrut_stable_aslr, ...
                        'cmip5_cf_aslr',cmip_counterfactual_aslr, ...
                        'budget',budget.central_aslr));
save('./data/source_data/fig1_data.mat','fig1dat')

%% Function Library

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