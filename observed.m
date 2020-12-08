% This code calculates the NY Battery and global observed trends in SLR
% directly from NOAA data (NY) and from Dangendorf et al. 2019 (global).
%
% Last Updated: 11/23/2020
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
% If you have any questions or comments about this code, please contact
% Daniel Gilford at dgilford@climatecentral.org
% ---
%

%% Setup
clear
close all

%% Load Data from NOAA

% NY Battery Tide Gauge (no. 8518750) data from NOAA's online archive:
% https://tidesandcurrents.noaa.gov/sltrends/sltrends_station.shtml?id=8518750
% Data accessed 29 July 2020

% define path to NOAA data file
dat_path='./data/observed/8518750_meantrend.csv';
% load the NOAA data as a table, convert to an array
T=readtable(dat_path,'HeaderLines',1);
dat=table2array(T(:,1:7));

%% Organize NOAA Data for Analysis

% get the time grid
nyear=dat(end,1)-dat(1,1)+1;
yeargrid=(1:1:nyear)'+dat(1,1)-1;
nmon=12;

% convert the data from month-over-month to yr/mon format
% third column is the Monthly Mean Sea Level
dat_ym=zeros(nyear,nmon);
for y=1:nyear
    for m=1:nmon
        tind=(y-1)*12+m;
        % data ends in June of 2020, so fill remaining slots with missing
        if (y==nyear) & (m>6)
            dat_ym(y,m)=NaN;
        else
            dat_ym(y,m)=dat(tind,3);
        end
    end
end

% take the annual average
dat_annual=nanmean(dat_ym,2);

% calculate the linear trend over 1900-2012
yrinds=find(yeargrid>1899 & yeargrid<2013);
[noaa_median_slope,S] = polyfit(yeargrid(yrinds), dat_annual(yrinds), 1);
[noaa_median_trendline,~] = polyval(noaa_median_slope,yeargrid(yrinds),S);
% Integrate linear trend over 1900-2012 to calculate total SLR
nyear_trends=length(yrinds)-1;
noaa_median_SLR=(noaa_median_slope(1)*nyear_trends)*100;

% NOAA uncertainty (over 112 years) drawn from
% "Sea Level Variations of the United States 1854-2006, NOAA Technical
% Report NOS CO-OPS 5"
% See: https://tidesandcurrents.noaa.gov/sltrends/sltrends_station.shtml?id=8518750
trendstd_112yr=0.176119403/1000;

% define the z-scores for converting to confidence intervals from stdev.
Z95=1.96; Z90=1.65;

% calculate the slopes for the 5th and 95th percentiles
upper_slope=stdtoX_plus(trendstd_112yr,noaa_median_slope(1),Z90);
lower_slope=stdtoX_minus(trendstd_112yr,noaa_median_slope(1),Z90);

% Integrate linear trends over 1900-2012 to calculate total SLR
% and convert to cm
noaa_upper_SLR=upper_slope*nyear_trends*100;
noaa_lower_SLR=lower_slope*nyear_trends*100;

% NOTE that NY data includes GIA effects in SLR since 1900. These are
% removed in the study

% plot
figure(1)
    plot(yeargrid(yrinds),dat_annual(yrinds),'LineWidth',3)
    hold on
    plot(yeargrid(yrinds),noaa_median_trendline,'k--','LineWidth',3)
    title('Relative Sea Level, 8518750 The Battery, New York')
    legend('Annual Mean','Linear Trend','Location','Best')
    xlabel('Year')
    ylabel('Meters')
    axis('tight')
hold off

%% Calculate Dangendorf et al. (2019) trends

% Taken from:
%   Dangendorf, S., Hay, C., Calafat, F.M. et al. 
%   Persistent acceleration in global sea-level rise since the 1960s. 
%   Nat. Clim. Chang. 9, 705?710 (2019). https://doi.org/10.1038/s41558-019-0531-8

% Median and uncertainty in trends from GMSL valid over 1900-2015
% We assume that this trend holds over the slightly shorter period of
% 1900-2012:
d19_median=1.6;  % mm/yr
d19_std=0.4;     % mm/yr

% calculate the slopes for the 5th and 95th percentiles
d19_upper=stdtoX_plus(d19_std,d19_median,Z90);
d19_lower=stdtoX_minus(d19_std,d19_median,Z90);

% Integrate linear trends over 1900-2012 to calculate total SLR
% and convert to cm
d19_median_SLR=d19_median*nyear_trends/10;
d19_upper_SLR=d19_upper*nyear_trends/10;
d19_lower_SLR=d19_lower*nyear_trends/10;

%% Aggregate and Print to Screen
ny_SLR_in_2012=[noaa_median_SLR,noaa_lower_SLR,noaa_upper_SLR] % includes GIA
global_SLR_in_2012=[d19_median_SLR,d19_lower_SLR,d19_upper_SLR]

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