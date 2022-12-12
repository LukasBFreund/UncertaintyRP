%%=========================================================================
% Main replication code for 
% "Unexpected Effects: Uncertainty, Unemployment, and Inflation"
% Freund, Lee, and Rendahl (2022)
%
% Baseline model
%
% Run on Matlab R2020b, Dynare 4.4.3
% Written: January 2019, last updated: December 2022
% For any questions please email
% lukas.beat.freund@gmail.com / hanbaeklee1@gmail.com
%%=========================================================================
 
%% Housekeeping
%--------------------------------------------------------------------------
clear; close all; clc; TimeStart = tic;
addpath("./mainFLR_others");

%% Running baseline model
%---------------------------------------------------------------------------
mainFLR_baseline;                % baseline model

%% Running other specifications (for analysis)
%---------------------------------------------------------------------------
cd("./mainFLR_others/");
dynare dynareFLR_noRP;           % no risk premium channel
dynare dynareFLR_flexible;       % flexible price
dynare dynareFLR_flexible_noRP;  % flexible price + no risk premium channel
dynare dynareFLR_demand;         % with demand shock
dynare dynareFLR_phipi05;        % Taylor rule parameter = 5               
dynare dynareFLR_phipi1000;      % Taylor rule parameter = 1000              

delete *.jnl *.log *.asv
cd("..");

%% Plot Figures 
%---------------------------------------------------------------------------
mainFLR_plot;   
