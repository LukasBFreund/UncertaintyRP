%%=========================================================================
% Replication code for "Unexpected Effects: Uncertainty, Unemployment, and
% Inflation" by Freund, Lee, and Rendahl (2022)
% Run on Matlab R2019b, Dynare 4.4.3
% Lukas B. Freund + Hanbaek Lee (Calibration)
% Written: January 2019, last updated: December 2022
% For any questions please email me at 
% lukas.beat.freund@gmail.com / hanbaeklee1@gmail.com
%%=========================================================================
 
%% Housekeeping
%--------------------------------------------------------------------------
clear;
close all;
clc;
TimeStart = tic;

%% Parameters

% Risk aversion
vPar.xi = 0.5;

% Wage elasticity
vPar.TZElasticity = 19.1;

% Patience          
vPar.betta = 0.9981;

% Vacancy posting cost (hbar=0.7)
vPar.kappa = 0.23643;

% Habit parameters
vPar.CCHabit1 = 0.996; 
vPar.CCHabit2 = 0.090; 

% Alternating offer parameter
vPar.chi = 0.7900;

% Policy
vPar.phi_pi = 1.5; 
vPar.phi_r = 0;
vPar.phi_y = 0;

% Productivity process
vPar.rho_z = 0.983;           
vPar.sigma_zbar = 0.0080;    
vPar.rho_sigma_z = 0.913;       
vPar.sigma_sigma_z = ((0.0853*0.01)/vPar.sigma_zbar)*vPar.sigma_zbar;

% Interest rate process
vPar.rho_e_R  = 0; 
vPar.sigma_RBar = 0;
vPar.rho_sigma_R = 0; 
vPar.sigma_sigma_R = 0; 

% Preference process
vPar.rho_q =  0; 
vPar.sigma_qBar = 0; 
vPar.rho_sigma_q = 0; 
vPar.sigma_sigma_q = 0; 

% Save parameters 
save Parameters vPar;
save ./mainFLR_others/Parameters vPar;

%% Run dynare 
%--------------------------------------------------------------------------
% Call dynare - KEY STEP
dynare dynareFLR noclearall


%% Unconditional moment (from simulation)
%-------------------------------------------
%key allocations
%-------------------------------------------
% Risk free return (RFR)
% Risk premium (RP)
% Sharpe ratio (=mean(RP)/std(RP))
% Price-dividend ratio (J/Dividend)

%PDratio
PDratio = J./Dividend;
%winsorize
PDratio = (PDratio<prctile(PDratio,1))*prctile(PDratio,1) + ...
          (PDratio>prctile(PDratio,99))*prctile(PDratio,99) + ...
          (PDratio>=prctile(PDratio,1)).*(PDratio<=prctile(PDratio,99)).*PDratio;

%excess return
R = ((1-delta)*J(2:end))./(J(1:(end-1))-(x(1:(end-1)).*z(1:(end-1))-w(1:(end-1))));
R = [R;R(1)];
excR = (R-RFR);
%winsorize
excR    = (excR<prctile(excR,1))*prctile(excR,1) + ...
          (excR>prctile(excR,99))*prctile(excR,99) + ...
          (excR>=prctile(excR,1)).*(excR<=prctile(excR,99)).*excR;
      
%-------------------------------------------
%quarterly ver.
%-------------------------------------------
%size of table
mTableQ1 = zeros(8,2);

%quarterly aggregation
numGroup = floor(length(u)/3);
group = [kron(1:numGroup,[1,1,1]),numGroup+1];
uQ = accumarray(group',u)/3;
yQ = accumarray(group',y);
RFR = accumarray(group',RFR-1);
excR = accumarray(group',excR);
[~,ufiltered] = hpfilter(log(uQ),1600);
[~,yfiltered] = hpfilter(log(yQ),1600);

dep     = (RFR)*100;
indep   = (RFR)*100;
coeff1  = regress(dep(2:end),indep(1:(end-1)));
dep     = excR*100;
indep   = excR*100;
coeff2  = regress(dep(2:end),indep(1:(end-1)));

mTableQ1(1,1:2) = [0.24,mean((RFR)*100)];         % mean risk-free return (%p.q) 
mTableQ1(2,1:2) = [0.40,std((RFR)*100)];          % volatility of risk-free return (%p.q) 
mTableQ1(3,1:2) = [0.88,coeff1(1)];               % AR(1) coefficient of risk-free return (p.q)
mTableQ1(4,1:2) = [1.47,mean(excR*100)];          % mean excess return (%p.q) 
mTableQ1(5,1:2) = [8.46,std(excR*100)];           % volatility of excess return (%p.q) 
mTableQ1(6,1:2) = [0.08,coeff2(1)];               % AR(1) coefficient of excess return (p.q)
mTableQ1(7,1:2) = [12.50,std(ufiltered)*100];     % Unemployment volatility
mTableQ1(8,1:2) = [2.06,std(yfiltered)*100];      % Output volatility

%-------------------------------------------
%report
%-------------------------------------------
fprintf(' \n');
fprintf(' \n');
fprintf(' \n');
fprintf('==============================\n');
fprintf('Calibration results (Table 2) \n');
fprintf('==============================\n');
fprintf(' \n');
disp(mTableQ1(1:6,:));
fprintf('------------------------------\n');
fprintf(' \n');
disp(mTableQ1(7:8,:));
fprintf('==============================\n');

%% Done!
%----------------------------------------------------------------------------
TimeEnd = toc(TimeStart);
fprintf(' \n');
disp(['Total run time was ',num2str(TimeEnd),' seconds']);

%% Clean
%----------------------------------------------------------------------------
delete *.jnl *.log *.asv dynareFLR_static.m dynareFLR_set_auxiliary_variables.m dynareFLR_dynamic.m dynareFLR_results.mat dynareFLR.m