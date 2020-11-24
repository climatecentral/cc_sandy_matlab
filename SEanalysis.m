% This code ...
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

% Naming conventions for data
% LinRate==Cooling Counterfactual, Mean==Stable Counterfactual
% Mar==Marcott, Mn==Mann
% HadCRUT4=Historical SL change
% Natural==CF=Counterfactual scenario for CMIP5

% define the quantile levels we want to calculate
q=[.05 .5 .95];

% define the year to look at statistics for (end of analysis, 1900-2012)
yrwant=2012;

%% Load and Organize HadCRUT4 SE analysis

% define HadCRUT4 filenames and paths
datpath='./data/SEanalysis/hadcrut/';
fn1={'Mar_LinRate.mat','Mar_Mean.mat','Mn_LinRate.mat','Mn_Mean.mat'};
fn2={'Mar_HadCRUT.mat','Mn_HadCRUT.mat'};

% pull in the time, parameters, and indices, and find the index we need
[dat_init]=load_SE_samples(fullfile(datpath,fn1{1}),'hadcrut');
% get the time grid
time=dat_init.time;
yrind=find(time==yrwant);
% parameter and temperature set indices
nP=dat_init.metadata.nP; % number of parameter sets for the semi-empirical model
Pset=dat_init.metadata.P_index;
nT=dat_init.metadata.nTemps_fixedP; % number of unique temperature histories the model is run with
clear dat_init % we no longer need this particular data structure
% define the array to remap between HadCRUT4 and counterfactual scenarios
refidx=[1 1 2 2];

% load the HadCRUT4 Historical data from each scenario (Marcott, Mann)
for qi=1:length(fn2)
    % get the data structure
    dati=load_SE_samples(fullfile(datpath,fn2{qi}),'hadcrut');
    % define the sea level rise relative to 1900
    sl2(:,:,qi)=dati.sl;
    clear dati % we no longer need this particular data structure
    % calculate the quantiles
    qsl2(:,:,qi)=quantile(sl2(:,:,qi),q);
end

% pool the Marcott+Mann samples to create the Historical summary statistics 
sl2_pooled=vertcat(sl2(:,:,1),sl2(:,:,2)); % pool the sea levels together
qsl2_pooled=quantile(sl2_pooled(:,:),q); % find the quantiles of the pooled samples

% load the HadCRUT4 Counterfactual Scenario data (stable,cooling/Marcott, Mann)
for qi=1:length(fn1)
    % get the data structure
    dati=load_SE_samples(fullfile(datpath,fn1{qi}),'hadcrut');
    % define the sea level rise relative to 1900
    sl1(:,:,qi)=dati.sl;
    clear dati % we no longer need this particular data structure
    % calculate the quantiles
    qsl1(:,:,qi)=quantile(sl1(:,:,qi),q);
    % calculate the differences between the historical and counterfactual
    % (this is the ASLR)
    sldiff(:,:,qi)=sl2(:,:,refidx(qi))-sl1(:,:,qi); % difference the historical and counterfactuals
    qsl1diff(:,:,qi)=quantile(sldiff(:,:,qi),q); % get the quantiles of the differences
end
% pool the samples to create the summary statistics
sl1_pooled(:,:,1)=vertcat(sl1(:,:,1),sl1(:,:,3)); % pool the sea levels together (cooling)
sl1_pooled(:,:,2)=vertcat(sl1(:,:,2),sl1(:,:,4)); % pool the sea levels together (stable)
qsl1_pooled=quantile(sl1_pooled(:,:,:),q); % find the quantiles of the pooled values

% take differences across samples, fixing temperature history and
% parameters
sldiff_pooled(:,:,1)=vertcat(sldiff(:,:,1),sldiff(:,:,3)); % pool the sea levels together (cooling)
sldiff_pooled(:,:,2)=vertcat(sldiff(:,:,2),sldiff(:,:,4)); % pool the sea levels together (stable)
qsl1diff_pooled=quantile(sldiff_pooled(:,:,:),q); % find the quantiles of the pooled differences

%% CMIP5

% define CMIP5 filepath
cmip5_datpath='./data/SEanalysis/cmip5/';

% now get ready to match identifiers to CMIP5 models
fid=fopen(fullfile(cmip5_datpath,'historical_CMIP35_index.dat')); % open the index file
C=textscan(fid,'%u %s %s %s','MultipleDelimsAsOne',1,'HeaderLines',12,'Delimiter',' '); % load the data
fclose(fid); % close the index file
% store the loaded data into a structure for use
CMIPHist.id=C{1};
CMIPHist.model=C{2};
CMIPHist.replicate=C{3};
CMIPHist.expt=C{4};
% repeat steps above with Natural forcing files
fid=fopen(fullfile(cmip5_datpath,'historicalNat_CMIP35_index.dat'));
C=textscan(fid,'%u %s %s %s','MultipleDelimsAsOne',1,'HeaderLines',12,'Delimiter',' ');
fclose(fid);
CMIPNat.id=C{1};
CMIPNat.model=C{2};
CMIPNat.replicate=C{3};
CMIPNat.expt=C{4};

% load samples from CMIP5 files
% counterfactual data
files3 = dir(fullfile(cmip5_datpath,'CF*.mat')); % index the filepath
Nsel=1000; Nt=length(time); % pull all 1000 samples that were created
sl3Mar=ones(Nsel,Nt,length(files3))*NaN; % init arrays
sl3Mn=ones(Nsel,Nt,length(files3))*NaN;
% loop over all files
for qqq=1:length(files3)
    % load the individual file in
    dofile=fullfile(cmip5_datpath,files3(qqq).name);
    [sl3,dotime,~,~,indices3,sp]=ImportKlausSESamples_ext(dofile,Nsel);
    [inttime,doi,doj]=intersect(dotime{1},time);
    colCMIP5(qqq)=sp.column_CMIP5;
    % store the sea levels from the CMIP5 runs
    sl3Mar(:,doj,qqq)=sl3{1}(:,doi);
    sl3Mn(:,doj,qqq)=sl3{2}(:,doi);
    % calculate the quantiles
    qsl3Mar(:,:,qqq)=quantile(sl3Mar(:,:,qqq),q);
    qsl3Mn(:,:,qqq)=quantile(sl3Mar(:,:,qqq),q);
    % get the indices
    sub=find(CMIPNat.id==colCMIP5(qqq));
    CFmodel{qqq}=CMIPNat.model{sub};
    CFreplicate{qqq}=CMIPNat.replicate{sub};
end
% Repeat above for historical data
files4 = dir(fullfile(cmip5_datpath,'CMIP5-hist*.mat'));
Nt=length(time);
sl4Mar=ones(Nsel,Nt,length(files4))*NaN;
sl4Mn=ones(Nsel,Nt,length(files4))*NaN;
for qqq=1:length(files4)
    dofile=fullfile(cmip5_datpath,files4(qqq).name);
    [sl4,dotime,~,~,indices4,sp]=ImportKlausSESamples_ext(dofile,indices3);
    [inttime,doi,doj]=intersect(dotime{1},time);
    colCMIP5_4(qqq)=sp.column_CMIP5;
    sl4Mar(:,doj,qqq)=sl4{1}(:,doi);
    sl4Mn(:,doj,qqq)=sl4{2}(:,doi);
    qsl4Mar(:,:,qqq)=quantile(sl4Mar(:,:,qqq),q);
    qsl4Mn(:,:,qqq)=quantile(sl4Mn(:,:,qqq),q);

    sub=find(CMIPHist.id==colCMIP5_4(qqq));
    Histmodel{qqq}=CMIPHist.model{sub};
    Histreplicate{qqq}=CMIPHist.replicate{sub};
end

% now match up Nat and Hist runs, and take differences
pooledmodel={CFmodel{:},Histmodel{:}} % print model names to screen
[u,ui,uj]=unique(pooledmodel);
CFmodelid=uj(1:length(CFmodel));
Histmodelid=uj((length(CFmodel)+1):end);
modelnames=u;
modelids=1:length(modelnames);

% take the differences, which define the CMIP5 global ASLR
count=1;
for www=1:length(modelids)
    disp(modelnames{www});
    subCF=find(CFmodelid==modelids(www));
    subHist=find(Histmodelid==modelids(www));
    if (length(subCF)>0)&&(length(subHist)>0)
        [doreplicates,dori,dorj]=intersect(CFreplicate(subCF),Histreplicate(subHist));
        subCF=subCF(dori); subHist=subHist(dorj);
        for qqq=1:length(subCF)
            % calculate and store the differences
            CMIPsldiffrep(count)=CFreplicate(subCF(qqq));
            CMIPsldiffmodel(count)=CFmodel(subCF(qqq));
            CMIPsldiffCFrun(count)=subCF(qqq);
            CMIPsldiffHistrun(count)=subHist(qqq);
            CMIPsldiffMar(:,:,count)=-sl3Mar(:,:,subCF(qqq))+sl4Mar(:,:,subHist(qqq));
            CMIPsldiffMn(:,:,count)=-sl3Mn(:,:,subCF(qqq))+sl4Mn(:,:,subHist(qqq));
            qCMIPsldiffMar(:,:,count)=quantile(CMIPsldiffMar(:,:,count),q);
            qCMIPsldiffMn(:,:,count)=quantile(CMIPsldiffMn(:,:,count),q);
            
            % store sea levels in addition to a difference
            cmipsl_cf_Mn(:,:,count)=sl3Mn(:,:,subCF(qqq));
            cmipsl_cf_Mar(:,:,count)=sl3Mar(:,:,subCF(qqq));
            cmipsl_hist_Mn(:,:,count)=sl4Mn(:,:,subHist(qqq));
            cmipsl_hist_Mar(:,:,count)=sl4Mar(:,:,subHist(qqq));
            
            count=count+1;
        end
        
    end
end

% with everything above, we could get a little lost. Let's reorganize the
% CMIP5 data before proceeding
nmodels=size(cmipsl_cf_Mar,3);
for m=1:nmodels
    % organize the data
    slcmip(:,:,m,1)=cmipsl_hist_Mar(:,:,m);
    slcmip(:,:,m,2)=cmipsl_cf_Mar(:,:,m);
    slcmip(:,:,m,3)=cmipsl_hist_Mn(:,:,m);
    slcmip(:,:,m,4)=cmipsl_cf_Mn(:,:,m);
    qcmip=quantile(slcmip,q); % find the quantiles of each individual model
    
    % take the differences along each models samples (ASLR)
    sldiffcmip(:,:,m,1)=slcmip(:,:,m,1)-slcmip(:,:,m,2); % difference the historical and counterfactuals (marcott)
    sldiffcmip(:,:,m,2)=slcmip(:,:,m,3)-slcmip(:,:,m,4); % difference the historical and counterfactuals (mann)
    qsldiffcmip=quantile(sldiffcmip,q); % find the quantiles of each individual model difference
end
% define the model names
model_names=CMIPsldiffmodel;
replicate_names=CMIPsldiffrep;

% drop out the GFDL model, which had faulty data (cc: Klaus Bittermann)
DROP=1;
if DROP==1
    dropname='GFDL-CM3'; dropind=find(contains(model_names,dropname));
    slcmip(:,:,dropind,:)=[]; qcmip(:,:,dropind,:)=[];
    sldiffcmip(:,:,dropind,:)=[]; qsldiffcmip(:,:,dropind,:)=[];
end

% pool the data across the models (reformat the data)
slcmip_all=nan(size(slcmip,1)*size(slcmip,3),size(slcmip,2),size(slcmip,4));
sldiffcmip_all=nan(size(slcmip,1)*size(slcmip,3),size(slcmip,2),2);
for r=1:size(slcmip,1)
    for m=1:size(slcmip,3)
        tind=(r-1)*m+m;
        slcmip_all(tind,:,:)=slcmip(r,:,m,:);
        sldiffcmip_all(tind,:,:)=sldiffcmip(r,:,m,:);
    end
end
% take quantiles across the pool of runs
qslcmip_all=quantile(slcmip_all,q);
qsldiffcmip_all=quantile(sldiffcmip_all,q);

% pool across scenarios (summary) 
slcmip_summary=nan(size(slcmip_all,1)*2,size(slcmip_all,2),2);
sldiffcmip_summary=nan(size(sldiffcmip_all,1)*size(sldiffcmip_all,3),size(slcmip_all,2));
for r2=1:size(slcmip_all,1)
    for s=1:2:3
        tind2=(r2-1)*ceil(mod(s/2,2))+ceil(mod(s/2,2));
        slcmip_summary(tind2,:,1)=slcmip_all(r2,:,s);
        slcmip_summary(tind2,:,2)=slcmip_all(r2,:,s+1);
        
        sind2=floor(mod(s/2,2))+1;
        sldiffcmip_summary(tind2,:)=sldiffcmip_all(r2,:,sind2);
    end
end

% % take quantiles across the full pool of runs and scenarios
qslcmip_summary=quantile(slcmip_summary,q);
qsldiffcmip_summary=quantile(sldiffcmip_summary,q);


%% Create structures for the sampled data (nested structure following...)

% With the full data load/analysis complete to calculate ASLR, we can store
% the results in a MATLAB structure for use in other codes and archiving

% define the metadata for the structures
metadata_semi=struct('quantiles',q,'final_year',yrwant,'Nparams_semi',nP,'NTemps_semi',nT,'Parameter_index',Pset);
metadata_cmip5=struct('quantiles',q,'final_year',yrwant,'Nmodels',size(slcmip,3),'Nparams_cmip5',size(slcmip,1), ...
                      'model_names',model_names,'replicate_names',replicate_names,'Parameter_index',Pset);
% define the information about the variables (names and units)
varinfo=struct('q',{'Quantiles'},'qdiff',{'Differences of Quantiles'}, ...
                'sl',{'Sea level relative to 1900 (cm)'}, ...
                'sldiff',{'Sea-level Differences between Historical and Counterfactual (cm)'}, ...
                'slq2012',{'Quantiles in yrwant (2012)'});
varinfo_cmip5=struct('q',{'Quantiles'},'qdiff',{'Differences of Quantiles'}, ...
                'sl',{'Sea level relative to 1900 (cm)'}, ...
                'sldiff',{'Sea-level Differences between Historical and Counterfactual (cm)'}, ...
                'slq2012',{'Quantiles in yrwant_cmip5 (2012)'});

% --------dataset--------
% hadcrut.
% cmip5.
% 
% --------scenarios--------
% .historical.
% .cf. [[cmip5 only]]
% .cooling. [[hadcrut only]]
% .stable. [[hadcrut only]]
% 
% --------columns--------
% .summary.
% .marcott.
% .mann.
% 
% --------results/variables--------
% .q=
% .qdiff=
% .sl=
% .sldiff=
% .slq2012=
%
% [[cmip5 only]]
% .qmodels=
% .qdiffmodels=
% .slmodels= 
% .sldiffmodels=

% create the hadcrut structure
hadcrut=struct( 'metadata', metadata_semi, ...
                'varinfo', varinfo, ...
                'historical', ...
                    struct('summary', ...
                        struct('q', squeeze(qsl2_pooled), ...
                               'qdiff', [NaN], ...
                               'sl', squeeze(sl2_pooled), ...
                               'sldiff', [NaN], ...
                               'slq2012', squeeze(qsl2_pooled(:,yrind))), ...
                     'marcott', ...
                        struct('q', squeeze(qsl2(:,:,1)), ...
                               'qdiff', [NaN], ...
                               'sl', squeeze(sl2(:,:,1)), ...
                               'sldiff', [NaN], ...
                               'slq2012', squeeze(qsl2(:,yrind,1))), ...
                    'mann', ...
                        struct('q', squeeze(qsl2(:,:,2)), ...
                               'qdiff', [NaN], ...
                               'sl', squeeze(sl2(:,:,2)), ...
                               'sldiff', [NaN], ...
                               'slq2012', squeeze(qsl2(:,yrind,2)))), ...
                'cooling', ... 
                    struct('summary', ...
                        struct('q', squeeze(qsl1_pooled(:,:,1)), ...
                               'qdiff', squeeze(qsl1diff_pooled(:,:,1)), ...
                               'sl', squeeze(sl1_pooled(:,:,1)), ...
                               'sldiff', squeeze(sldiff_pooled(:,:,1)), ...
                               'slq2012', squeeze(qsl1_pooled(:,yrind,1))), ...
                     'marcott', ...
                        struct('q', squeeze(qsl1(:,:,1)), ...
                               'qdiff', squeeze(qsl1diff(:,:,1)), ...
                               'sl', squeeze(sl1(:,:,1)), ...
                               'sldiff', squeeze(sldiff(:,:,1)), ...
                               'slq2012', squeeze(qsl1(:,yrind,1))), ...
                    'mann', ...
                        struct('q', squeeze(qsl1(:,:,3)), ...
                               'qdiff', squeeze(qsl1diff(:,:,3)), ...
                               'sl', squeeze(sl1(:,:,3)), ...
                               'slqdiff', squeeze(sldiff(:,:,3)), ...
                               'slq2012', squeeze(qsl1(:,yrind,3)))), ...
                'stable', ...
                    struct('summary', ...
                        struct('q', squeeze(qsl1_pooled(:,:,2)), ...
                               'qdiff', squeeze(qsl1diff_pooled(:,:,2)), ...
                               'sl', squeeze(sl1_pooled(:,:,2)), ...
                               'sldiff', squeeze(sldiff_pooled(:,:,2)), ...
                               'slq2012',  squeeze(qsl1_pooled(:,yrind,2))), ...
                     'marcott', ...
                        struct('q', squeeze(qsl1(:,:,2)), ...
                               'qdiff', squeeze(qsl1diff(:,:,2)), ...
                               'sl', squeeze(sl1(:,:,2)), ...
                               'slqdiff', squeeze(sldiff(:,:,2)), ...
                               'sl2012', squeeze(qsl1(:,yrind,2))), ...
                     'mann', ...
                        struct('q', squeeze(qsl1(:,:,4)), ...
                               'qdiff', squeeze(qsl1diff(:,:,4)), ...
                               'sl', squeeze(sl1(:,:,4)), ...
                               'slqdiff', squeeze(sldiff(:,:,4)), ...
                               'slq2012', squeeze(qsl1(:,yrind,4)))) ...
    )


% create the cmip5 structure
% recall:
% (scenarios)       1=Marcott Hist, 2=Marcott CF, 3=Mann Hist, 4=Mann CF
% (history)         1=Marcott, 2=Mann
% (summary/history) 1=Hist, 2=CF
% (arrays to structure)
%   q       --> qslcmip_all(:,:,scenarios/history)
%   qdiff   --> qsldiffcmip_all(:,:,scenarios)
%   sl      --> slcmip_all(:,:,scenarios/history)
%   sldiff  --> sldiffcmip_all(:,:,scenarios)
%   slq2012 --> qslcmip_all(:,yrind,scenarios/history)
%   qmodels      --> qcmip(:,:,:,scenarios/history)
%   qdiffmodels  --> qsldiffcmip(:,:,:,scenarios)
%   slmodels     --> slcmip(:,:,:,scenarios/history)
%   sldiffmodels --> sldiffcmip(:,:,scenarios)
% SUMMARY
%   q       --> qslcmip_summary(:,:,history)
%   qdiff   --> qsldiffcmip_summary(:,:)
%   sl      --> slcmip_summary(:,:,history)
%   sldiff  --> sldiffcmip_summary(:,:)
%   slq2012 --> qslcmip_summary(:,yrind,history)
%   qmodels      --> [NaN]
%   qdiffmodels  --> [NaN]
%   slmodels     --> [NaN]
%   sldiffmodels --> [NaN]

cmip5=struct( 'metadata', metadata_cmip5, ...
                'varinfo', varinfo_cmip5, ...
                'historical', ...
                    struct('summary', ...
                        struct('q', qslcmip_summary(:,:,1), ...
                               'qdiff', [NaN], ...
                               'sl', slcmip_summary(:,:,1), ...
                               'sldiff', [NaN], ...
                               'slq2012', qslcmip_summary(:,yrind,1), ...
                               'qmodels', [NaN], ...
                               'qdiffmodels', [NaN], ...
                               'slmodels', [NaN], ...
                               'sldiffmodels', [NaN]), ...
                     'marcott', ...
                        struct('q', qslcmip_all(:,:,1), ...
                               'qdiff', [NaN], ...
                               'sl', slcmip_all(:,:,1), ...
                               'sldiff', [NaN], ...
                               'slq2012', qslcmip_all(:,yrind,1), ...
                               'qmodels', qcmip(:,:,:,1), ...
                               'qdiffmodels', [NaN], ...
                               'slmodels', slcmip(:,:,:,1), ...
                               'sldiffmodels', [NaN]), ...
                    'mann', ...
                        struct('q', qslcmip_all(:,:,3), ...
                               'qdiff', [NaN], ...
                               'sl', slcmip_all(:,:,3), ...
                               'sldiff', [NaN], ...
                               'slq2012', qslcmip_all(:,yrind,3), ...
                               'qmodels', qcmip(:,:,:,3), ...
                               'qdiffmodels', [NaN], ...
                               'slmodels', slcmip(:,:,:,3), ...
                               'sldiffmodels', [NaN])), ...
                'cf', ... 
                    struct('summary', ...
                        struct('q', qslcmip_summary(:,:,2), ...
                               'qdiff', qsldiffcmip_summary, ...
                               'sl', slcmip_summary(:,:,2), ...
                               'sldiff', sldiffcmip_summary, ...
                               'slq2012', qslcmip_summary(:,yrind,2), ...
                               'qmodels', [NaN], ...
                               'qdiffmodels', [NaN], ...
                               'slmodels', [NaN], ...
                               'sldiffmodels', [NaN]), ...
                     'marcott', ...
                        struct('q', qslcmip_all(:,:,2), ...
                               'qdiff', qsldiffcmip_all(:,:,1), ...
                               'sl', slcmip_all(:,:,2), ...
                               'sldiff', sldiffcmip_all(:,:,1), ...
                               'slq2012', qslcmip_all(:,yrind,2), ...
                               'qmodels', qcmip(:,:,:,2), ...
                               'qdiffmodels', qsldiffcmip(:,:,:,1), ...
                               'slmodels', slcmip(:,:,:,2), ...
                               'sldiffmodels', sldiffcmip(:,:,1)), ...
                    'mann', ...
                        struct('q', qslcmip_all(:,:,4), ...
                               'qdiff', qsldiffcmip_all(:,:,2), ...
                               'sl', slcmip_all(:,:,4), ...
                               'sldiff', sldiffcmip_all(:,:,2), ...
                               'slq2012', qslcmip_all(:,yrind,4), ...
                               'qmodels', qcmip(:,:,:,4), ...
                               'qdiffmodels', qsldiffcmip(:,:,:,2), ...
                               'slmodels', slcmip(:,:,:,4), ...
                               'sldiffmodels', sldiffcmip(:,:,2))) ...
    )

% save the results
save('./data/fig1/SEanalysis.mat','cmip5','hadcrut');

%% Create a table with the data

% generate the column names
table1_columns={'Scenario','Summary50','Summary5', 'Summary95','Mann50','Mann5', 'Mann95','Marcott50','Marcott5', 'Marcott95'};

% CMIP5 simulations
% generate the rows
CMIP5_based_rows={'Historical','Counterfactual'}';
CMIP5_based_simulations=table(CMIP5_based_rows, ...
    [round(cmip5.historical.summary.q(2,yrind),1); round(cmip5.cf.summary.q(2,yrind),1)], ...
    [round(cmip5.historical.summary.q(1,yrind),1); round(cmip5.cf.summary.q(1,yrind),1)], ...
    [round(cmip5.historical.summary.q(3,yrind),1); round(cmip5.cf.summary.q(3,yrind),1)], ...
    [round(cmip5.historical.mann.q(2,yrind),1); round(cmip5.cf.mann.q(2,yrind),1)], ...
    [round(cmip5.historical.mann.q(1,yrind),1); round(cmip5.cf.mann.q(1,yrind),1)], ...
    [round(cmip5.historical.mann.q(3,yrind),1); round(cmip5.cf.mann.q(3,yrind),1)], ...
    [round(cmip5.historical.marcott.q(2,yrind),1); round(cmip5.cf.marcott.q(2,yrind),1)], ...
    [round(cmip5.historical.marcott.q(1,yrind),1); round(cmip5.cf.marcott.q(1,yrind),1)], ...
    [round(cmip5.historical.marcott.q(3,yrind),1); round(cmip5.cf.marcott.q(3,yrind),1)], ...
    'VariableNames',table1_columns)

% HadCRUT4 simulations
% generate the rows
T_based_rows={'Historical','Stable','Cooling'}';
T_based_simulations=table(T_based_rows, ...
    [round(hadcrut.historical.summary.q(2,yrind),1); round(hadcrut.stable.summary.q(2,yrind),1); round(hadcrut.cooling.summary.q(2,yrind),1)], ...
    [round(hadcrut.historical.summary.q(1,yrind),1); round(hadcrut.stable.summary.q(1,yrind),1); round(hadcrut.cooling.summary.q(1,yrind),1)], ...
    [round(hadcrut.historical.summary.q(3,yrind),1); round(hadcrut.stable.summary.q(3,yrind),1); round(hadcrut.cooling.summary.q(3,yrind),1)], ...
    [round(hadcrut.historical.mann.q(2,yrind),1); round(hadcrut.stable.mann.q(2,yrind),1); round(hadcrut.cooling.mann.q(2,yrind),1)], ...
    [round(hadcrut.historical.mann.q(1,yrind),1); round(hadcrut.stable.mann.q(1,yrind),1); round(hadcrut.cooling.mann.q(1,yrind),1)], ...
    [round(hadcrut.historical.mann.q(3,yrind),1); round(hadcrut.stable.mann.q(3,yrind),1); round(hadcrut.cooling.mann.q(3,yrind),1)], ...
    [round(hadcrut.historical.marcott.q(2,yrind),1); round(hadcrut.stable.marcott.q(2,yrind),1); round(hadcrut.cooling.marcott.q(2,yrind),1)], ...
    [round(hadcrut.historical.marcott.q(1,yrind),1); round(hadcrut.stable.marcott.q(1,yrind),1); round(hadcrut.cooling.marcott.q(1,yrind),1)], ...
    [round(hadcrut.historical.marcott.q(3,yrind),1); round(hadcrut.stable.marcott.q(3,yrind),1); round(hadcrut.cooling.marcott.q(3,yrind),1)], ...
    'VariableNames',table1_columns)

% Differences from historical (ASLR)
Diff_rows={'CMIP5','Stable','Cooling'}';
Diff_from_Hist_ASLR=table(Diff_rows, ...
    [round(cmip5.cf.summary.qdiff(2,yrind),1); round(hadcrut.stable.summary.qdiff(2,yrind),1); round(hadcrut.cooling.summary.qdiff(2,yrind),1)], ...
    [round(cmip5.cf.summary.qdiff(1,yrind),1); round(hadcrut.stable.summary.qdiff(1,yrind),1); round(hadcrut.cooling.summary.qdiff(1,yrind),1)], ...
    [round(cmip5.cf.summary.qdiff(3,yrind),1); round(hadcrut.stable.summary.qdiff(3,yrind),1); round(hadcrut.cooling.summary.qdiff(3,yrind),1)], ...
    [round(cmip5.cf.mann.qdiff(2,yrind),1); round(hadcrut.stable.mann.qdiff(2,yrind),1); round(hadcrut.cooling.mann.qdiff(2,yrind),1)], ...
    [round(cmip5.cf.mann.qdiff(1,yrind),1); round(hadcrut.stable.mann.qdiff(1,yrind),1); round(hadcrut.cooling.mann.qdiff(1,yrind),1)], ...
    [round(cmip5.cf.mann.qdiff(3,yrind),1); round(hadcrut.stable.mann.qdiff(3,yrind),1); round(hadcrut.cooling.mann.qdiff(3,yrind),1)], ...
    [round(cmip5.cf.marcott.qdiff(2,yrind),1); round(hadcrut.stable.marcott.qdiff(2,yrind),1); round(hadcrut.cooling.marcott.qdiff(2,yrind),1)], ...
    [round(cmip5.cf.marcott.qdiff(1,yrind),1); round(hadcrut.stable.marcott.qdiff(1,yrind),1); round(hadcrut.cooling.marcott.qdiff(1,yrind),1)], ...
    [round(cmip5.cf.marcott.qdiff(3,yrind),1); round(hadcrut.stable.marcott.qdiff(3,yrind),1); round(hadcrut.cooling.marcott.qdiff(3,yrind),1)], ...
    'VariableNames',table1_columns)

%% Function Library

function dat=load_SE_samples(file_in,type)
    % this code loads in the samples from the semi-empirical output files
    % from the semi-empirical model runs. Here we are only concerned with
    % the sea level output, and we have defined a fixed parameter set array
    % which is provided to the user
    % "type" is whether we are loading hadcrut or cmip5 files

    % load in the data file, turn off warninings for the load
    warning('off')
    fid=load(file_in);
    
    % find the number of samples, and total number of samples
    if strcmp(type,'hadcrut')==1
        settings=fid.S.settings;
        S=fid.S;
    elseif strcmp(type,'cmip5')==1
        settings=fid.S.Mar.settings;
        S=fid.S.Mar;
    else
        error('bad input type of file to load')
    end
    numsamps=settings.Tnum*settings.sample;
        metadata=struct('nsamps_total',numsamps, ...
            'nTemps_fixedP',settings.Tnum, ...
            'nP',settings.sample);

    % define the parameter index array (integers)
    p_index=1:1:settings.sample;
    metadata.P_index=p_index;

    % get the sea levels and temperatures
    dat.sl=S.sl;
    dat.T01=S.T01;
    dat.Tcf=S.Tcf;
    % grab the time index
    dat.time=S.time;
    dat.metadata=metadata;

    return
end

function [sl,t,temp,temp0,indices,specs]=ImportKlausSESamples_ext(fn,indices)
% Function to load in Klaus Bittermann's SES CMIP5 runs
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Jun 15 17:10:48 EDT 2016
y=load(fn);

    if isfield(y.S,'sl')

        defval('indices',1:size(y.S.sl,1));
        if length(indices)==1
            indices=round(linspace(2,size(y.S.sl,1),indices));
        end

        sl=y.S.sl(indices,:);
        t=y.S.time;
        temp=y.S.Tcf(indices,:);
        temp0=y.S.T01(indices,:);
        specs.fnames='';
    else
        fnames=fieldnames(y.S);
        specs.fnames=fnames;
        defval('indices',1:size(y.S.(fnames{1}).sl,1));

        if length(indices)==1
            indices=round(linspace(2,size(y.S.(fnames{1}).sl,1),indices));
        end
        for www=1:length(fnames)
            sl{www}=y.S.(fnames{www}).sl(indices,:);
            t{www}=y.S.(fnames{www}).time;
            temp{www}=y.S.(fnames{www}).Tcf(indices,:);
            temp0{www}=y.S.(fnames{www}).T01(indices,:);
        end
        fnames2=fieldnames(y.S.(fnames{1}));
        for ppp=1:length(fnames2)
            if length(y.S.(fnames{1}).(fnames2{ppp}))==1
                specs.(fnames2{ppp})=y.S.(fnames{1}).(fnames2{ppp});
            end
        end
    end
end

function defval(name,value)
% DEFVAL(name,value)
%
% Assigns a default value to the named variable
%
% INPUT:
% 
% name    A string, enclosed in single quotes, with a variable name
% value   The value, whatever it is, that you want the variable to have 
%
% OUTPUT:
%
%      None. The variables appear as if by magic into your workspace or
%      will be available inside your function.
%
% NOTE: 
%
% This won't work for an unassigned structure variable.
%
% Last modified by ebrevdo-at-alumni-princeton.edu, 05/28/2011
% Last modified by fjsimons-at-alum.mit.edu, 12/20/2012
% 
% It appears that defval('bla',functioncall) evaluates the function call
% regardless of whether or not 'bla' has been assigned.

    if ~ischar(name),
      error(sprintf(['The first argument of DEFVAL ',...
            'has to be a string with a variable name']));
    end

    % Always do it is our default here
    si=1;
    % If it exists... as a variable (say it, it makes it faster!)
    if evalin('caller',[ 'exist(''' name ''',''var'')']);
      % ... and it's empty, do it; but don't do it if it's non empty
      si=evalin('caller',[ 'isempty(' name ')']);
    end
    % Do it or not
    if si
      assignin('caller',name,value);
    end

end