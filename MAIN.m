%Main.m
% replication code of the final Assignment


clear all; close all; clc;

%% IMPORT THE DATA used in the analysis (Direct from the paper)
load ('DATA/CGV_quarterly_macro.mat');       % all vars are same length, but there are missing values. 
run DATA/GZ_spreads.m                   % Import data from the AER database of Gilchrist and Zakrajsek

% ----Quarterly measures, constructed using patents.csv file in QuarterlyXI.do; then run construct_patentindex.m:
load ('DATA/patentindex.mat');                 %lxiq_pc

%---- Load Industry DATA:
load ('DATA/xiq_industry.mat');           
load ('DATA/xiq_manufacturing.mat');      
load ('DATA/xiq_services.mat'); 
load ('DATA/ltfpf.mat')                 % unadjusted TFP
load ('DATA/EXT_SHOCKS.mat')
% ----------------------------------------------------------------------------------------
% Time Period:
%  All the data (with some NaN) are from 1947, but cut it, start from 1948.
begin_date      = 1948;  % 
end_date        = 2010;
% -------
samp  = find(ddate==1948):find(ddate==2010.75); % long sample
samp1 = find(ddate==1961):find(ddate==2010.75); % short sample
% samp2 = find(ddate==1973):find(ddate==2010.75); % sample (credit conditions)

% Choose for the proxy VAR:
instrument_var    = 12;
ext_shocks_ind    = 25:29;%
rep_proxy         = 1000000; % Number of draws in MC chain (takes time)

% Define Data:
DATA_MASTER(:,1)   = lvxo(samp);
DATA_MASTER(:,2)   = ly1(samp);                % Real log per capita output 
DATA_MASTER(:,3)   = lcnds(samp);              % Real log per capita consumption non durables and services  
DATA_MASTER(:,4)   = linv(samp);               % Real log per capita investment
DATA_MASTER(:,5)   = lh(samp);                 % Real log per capita hours
DATA_MASTER(:,6)   = lGDPdefl(samp);           % log GDP deflator
DATA_MASTER(:,7)   = ltfp_ua(samp);            % Log Fernald's TFP
DATA_MASTER(:,8)   = log(ics(samp));           % Log Consumer Confidence (Michigan Survey)
DATA_MASTER(:,9)   = infla(samp)./100;         % Inflation 
DATA_MASTER(:,10)  = ffr(samp)./100;           % Federal Funds rate
DATA_MASTER(:,11)  = lsp_pc(samp);             % Real log per capita stock price index
DATA_MASTER(:,12)  = lxiq_pc(samp);            % Real log per capita patent-based innovation index
DATA_MASTER(:,13)  = ltfpf(samp);              % Log Fernald's unadjusted TFP
DATA_MASTER(:,14)  = lxiq_industry_pc(samp,1); % Patent-based index: finance
DATA_MASTER(:,15)  = lxiq_industry_pc(samp,2); % Patent-based index: manufacturing
DATA_MASTER(:,16)  = lxiq_industry_pc(samp,3); % Patent-based index: mining
DATA_MASTER(:,17)  = lxiq_industry_pc(samp,4); % Patent-based index: services
DATA_MASTER(:,18)  = lxiq_industry_pc(samp,5); % Patent-based index: transportation
DATA_MASTER(:,19)  = lxiq_industry_pc(samp,6); % Patent-based index: wholesale
DATA_MASTER(:,20)  = GZ_spread(samp,1);         % Credit Conditions: GZ Credit spread (Gilchrist & Zakrajsek, AER 2012)
DATA_MASTER(:,21)  = GZ_spread(samp,5);         % Credit Conditions: EBP (Excess Bond Premium, Gilchris & Zakrajsek, AER 2012)
DATA_MASTER(:,22)  = GZ_spread(samp,2);         % Credit Conditions: BAA-AAA (from Gilchrist & Zakrajsek, AER 2012 database)
DATA_MASTER(:,23)  = lxiq_manuf_pc(samp,17);    % Patent-based index: Electronic And Other Electrical Equipment And Components, Except Computer Equipment
DATA_MASTER(:,24)  = lxiq_services_pc(samp,1);  % Patent-based index: Business Services
DATA_MASTER(:,25)  = EXT_SHOCKS(samp,1);        % News series about tax shocks from Leeper, Walker, and Yang (2012)
DATA_MASTER(:,26)  = EXT_SHOCKS(samp,2);        % Nominal present value of Ramey (2011) defense news variable divided by nominal gdp of previous quarter
DATA_MASTER(:,27)  = EXT_SHOCKS(samp,3);        % Hamilton (2003) net oil price increase (3 year)
DATA_MASTER(:,28)  = EXT_SHOCKS(samp,4);        % Romer and Romer (2004) monetary policy shocks constructed by Barakchian Crowe (2013) through 2006
DATA_MASTER(:,29)  = EXT_SHOCKS(samp,5);        % Mertens and Ravn (2012) unanticipated tax shock series
% Assign names to all the series in the dataset:
labels_vars =   [       '         VXO                  '; %1        
                        '      Output                  '; %2
                        '       Consumption            '; %3
                        '       Investment             '; %4
                        '        Hours                 '; %5
                        '      Price Level             '; %6			                                                 
                        '  Utilization-adjusted TFP    '; %7
                        '    Consumer Confidence       '; %8
                        '     Inflation                '; %9
                        ' Federal Funds Rate           '; %10
                        ' Stock Prices                 '; %11
                        ' Patent-Based Innovation Index'; %12
                        '       TFP  not adjusted      '; %13
                        ' Industry xiq -Finance        '; %14
                        'Industry xiq - Manuf          '; %15
                        'Industry xiq - Mining         '; %16
                        'Industry xiq - services       '; %17
                        'Industry xiq - transport      '; %18
                        'Industry xiq - wholesale      '; %19
                        '    GZ Credit Spread          '; %20
                        ' Excess Bond Premium          '; %21
                        ' BAA-AAA Spread               '; %22
                        '   Electric/Electronic        '; %23                      
                        ' Business Services            '; %24
                        ' News about tax               '; %25
                        ' News about govt. spending    '; %26
                        ' Oil price shock              '; %27
                        ' Monetary policy shock        '; %28
                        ' Tax shock                    '];%29 

nlags   = 4;                    % VAR(nlags) 
nconst  = 1;                   
prior   = 1;                                                                    
ndraws  = 1000;                  % number of draws
bound   = 0.16;                  % lower bound of confidence interval
                                % (e.g. for bound=0.16, CI is 16% - 84%)                                                              
nirf    = 21;                    % IRFs horizon                        

esty1 = 1961;     
estq1 = 1;        
esty2 = 2010;     
estq2 = 4;        

%% Replicate Figure 1 - KPSS index:
patentindex = lxiq_pc(samp1);          
%====
FIG1=figure('Color',[0.9412 0.9412 0.9412],'Position',[1 1 800-100 600-100],'Name','FIG1_CGV_RESTAT');   

figure(FIG1)
daux = [1961:0.25:2010.75];
patentindex= [lxiq_pc(5+(1961-1948)*4:end)]';
NBERbc(daux, patentindex, '-', 2, 'blue');
title ('Figure 1: Quarterly Patent-Based Innovation Index')
xlim ([1961 2010])
ylim ([-4 -1])
clear samp1



%% Replicate result Bayesian VAR
run Repl_Figs_BVAR.m 
% 
% 
%% Replicate result Bayesian Proxy VAR
run Repl_Figs_BPSVAR.m
%
%% Robustness Check: Local Projections
run Repl_Figs_LP.m 
%% Robustness Check:Bayesian Local Projections
run Repl_Figs_BLP.m 