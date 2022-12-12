%%=========================================================================
% 
% This file plots all the figures 
%
% Last updated: December 2022
%
%%=========================================================================

clear; close all; clc;

%--------------------------------------------------------------------------
%% Design choices
%--------------------------------------------------------------------------

% Choose whether to print and, if so, set target folder 
OptionPrint = 1;        
TargetPath ='./Output/Figures/';

% General design choices 
IRFPeriodsNum = 24;
optionGreycolor = 1;    
vLinestyle = {'-','--',':','-.'};
fontsizeDefault = 10;
fontsizeAxis = 10;
fontsizeAxisticks = 8;
fontsizeLegend = 8;
fonttype = 'times';
linewidthDefault = 2;
colorZeros = 'k';
styleZeros = ':';
fontweight = 'normal';

% Colors choices
if optionGreycolor == 0
    vColors = {[0,42,91]/255,[255 127 0]/255,[114,47,55]/255,[128,128,128]/255};
elseif optionGreycolor == 1
    vColors = {[0.2,0.2,0.2],[0.33 0.33 0.33],[0.46,0.46,0.46],[0.6,0.6,0.6]};
end

% Adjust some style options
set(groot, 'DefaultTextInterpreter', 'none');  
set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','none');


%--------------------------------------------------------------------------
%% Figure 1
%-------------------------------------------------------------------------

name = 'fig_1';

% Settings
NumModels = 2; 
vVariables = {'J','u','RFR','RP'}; 
vVNames = {'Firm value','Unemployment rate','Risk-free rate (ann.)','Risk premium (ann.)'};
NumSubplotV = 2;
NumSubplotH = 2;
OptionLegend = 1;
legend_labels{1} = 'Baseline';
legend_labels{2} = 'Without RP channel';
LegendPosition = 'northeast';

% Load data and, where applicable, transform
load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_flexible'));
mIRFProp_1= mIRFProp_zUncertainty_EMAS; 

uPos = strmatch('u',vNames,'exact');
PiiPos = strmatch('Pii',vNames,'exact');
RFRPos = strmatch('RFR',vNames,'exact');
RPPos = strmatch('RP',vNames,'exact');

mIRFProp_1(:,uPos,:) = mIRFProp_1(:,uPos,:)*vEMAS(uPos);
mIRFProp_1(:,PiiPos,:) = 12* mIRFProp_1(:,PiiPos,:)*vEMAS(PiiPos); 
mIRFProp_1(:,RFRPos,:) = 12*mIRFProp_1(:,RFRPos,:)*vEMAS(RFRPos);
mIRFProp_1(:,RPPos,:) = 100*mIRFProp_1(:,RPPos,:)*vEMAS(RPPos);

load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_flexible_noRP'));
mIRFProp_2= mIRFProp_zUncertainty_EMAS; 
mIRFProp_2(:,uPos,:) = mIRFProp_2(:,uPos,:)*vEMAS(uPos);
mIRFProp_2(:,PiiPos,:) = 12* mIRFProp_2(:,PiiPos,:)*vEMAS(PiiPos); 
mIRFProp_2(:,RFRPos,:) = 12*mIRFProp_2(:,RFRPos,:)*vEMAS(RFRPos);
mIRFProp_2(:,RPPos,:) = 100*mIRFProp_2(:,RPPos,:)*vEMAS(RPPos);

% Print
mainFLR_print_IRF;
close;


%--------------------------------------------------------------------------
%% Figure 2
%--------------------------------------------------------------------------

name = 'fig_2';
NumModels = 2; 
vVariables = {'J','u', 'c','x','Pii','RFR'};  
vVNames = {'Firm value','Unemployment rate','Consumption','Relative price','Inflation rate (ann.)','Risk-free rate (ann.)'};
NumSubplotV = 3;
NumSubplotH = 2;
OptionLegend = 1;
legend_labels{1} = 'Baseline';
legend_labels{2} = 'Without RP channel';
LegendPosition = 'southeast';

load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_baseline'));
mIRFProp_1= mIRFProp_zUncertainty_EMAS; 
uPos = strmatch('u',vNames,'exact');
PiiPos = strmatch('Pii',vNames,'exact');
RFRPos = strmatch('RFR',vNames,'exact');
RPPos = strmatch('RP',vNames,'exact');
mIRFProp_1(:,uPos,:) = mIRFProp_1(:,uPos,:)*vEMAS(uPos);
mIRFProp_1(:,PiiPos,:) = 12* mIRFProp_1(:,PiiPos,:)*vEMAS(PiiPos); %
mIRFProp_1(:,RFRPos,:) = 12*mIRFProp_1(:,RFRPos,:)*vEMAS(RFRPos);
mIRFProp_1(:,RPPos,:) = 100*mIRFProp_1(:,RPPos,:)*vEMAS(RPPos);

load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_noRP'));
mIRFProp_2= mIRFProp_zUncertainty_EMAS; 
mIRFProp_2(:,uPos,:) = mIRFProp_2(:,uPos,:)*vEMAS(uPos);
mIRFProp_2(:,PiiPos,:) = 12* mIRFProp_2(:,PiiPos,:)*vEMAS(PiiPos); 
mIRFProp_2(:,RFRPos,:) = 12*mIRFProp_2(:,RFRPos,:)*vEMAS(RFRPos);
mIRFProp_2(:,RPPos,:) = 100*mIRFProp_2(:,RPPos,:)*vEMAS(RPPos);

mainFLR_print_IRF;
close;
clear legend_labels;

%--------------------------------------------------------------------------
%% Figure 3
%--------------------------------------------------------------------------

name = 'fig_3';
mainFLR_plot_decomp;
close;

%--------------------------------------------------------------------------
%% Figure 4
%--------------------------------------------------------------------------

name = 'fig_4';
NumModels = 2; 
vVariables = {'Pii','RP'};
vVNames = {'Inflation rate (ann.)','Risk premium (ann.)'};
NumSubplotV = 1;
NumSubplotH = 2;
OptionLegend = 1;
legend_labels{1} = 'Uncertainty shock';
legend_labels{2} = 'Demand shock';
LegendPosition = 'southeast';

load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_demand'));
mIRFProp_1= mIRFProp_zUncertainty_EMAS; 
mIRFProp_2= mIRFProp_Andreasen_R_EMAS;  % standard GIRF for level shock 
 
uPos = strmatch('u',vNames,'exact');
PiiPos = strmatch('Pii',vNames,'exact');
RPPos = strmatch('RP',vNames,'exact');

mIRFProp_1(:,uPos,:) = mIRFProp_1(:,uPos,:)*vEMAS(uPos);
mIRFProp_2(:,uPos,:) = mIRFProp_2(:,uPos,:)*vEMAS(uPos);
mIRFProp_1(:,PiiPos,:) = 12* mIRFProp_1(:,PiiPos,:)*vEMAS(PiiPos); 
mIRFProp_2(:,PiiPos,:) = 12* mIRFProp_2(:,PiiPos,:)*vEMAS(PiiPos); 
mIRFProp_1(:,RPPos,:) = 100*mIRFProp_1(:,RPPos,:)*vEMAS(RPPos);
mIRFProp_2(:,RPPos,:) = 100*mIRFProp_2(:,RPPos,:)*vEMAS(RPPos);

% Linear scaling 
mIRFProp_2 = mIRFProp_2/(mIRFProp_2(1,uPos)/mIRFProp_1(1,uPos)); 

mainFLR_print_IRF;
close;

%--------------------------------------------------------------------------
%% Figure 5
%--------------------------------------------------------------------------

name = 'fig_5';
     
NumModels = 4; 
vVariables = {'u'};
vVNames = {'Unemployment rate'};
NumSubplotV = 1;
NumSubplotH = 1;
OptionLegend = 1;
legend_labels{1} = 'Baseline';
legend_labels{2} = '$\phi_{\pi} = 5$';
legend_labels{3} = '$\phi_{\pi} \rightarrow \infty$';
legend_labels{4} = 'Flex price';
LegendPosition = 'northeast';

load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_baseline')); 
mIRFProp_1= mIRFProp_zUncertainty_EMAS; 
uPos = strmatch('u',vNames,'exact');
mIRFProp_1(:,uPos,:) = mIRFProp_1(:,uPos,:)*vEMAS(uPos);

load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_phipi05')); 
mIRFProp_2= mIRFProp_zUncertainty_EMAS;  % use model-appropriate EMAS
mIRFProp_2(:,uPos,:) = mIRFProp_2(:,uPos,:)*vEMAS(uPos);

load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_phipi1000')); 
mIRFProp_3= mIRFProp_zUncertainty_EMAS; 
mIRFProp_3(:,uPos,:) = mIRFProp_3(:,uPos,:)*vEMAS(uPos);

load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_flexible')); 
mIRFProp_4= mIRFProp_zUncertainty_EMAS; 
mIRFProp_4(:,uPos,:) = mIRFProp_4(:,uPos,:)*vEMAS(uPos);

mainFLR_print_IRF;
close;
clear legend_labels;

close;






