function [ys,check] = dynareFLR_noRP_steadystate(ys,exo)

global M_ 
% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

%% Enter model equations here
    e_R = 0;
    z = zss;
    Pii = Piiss;
    h = hss; 
    u = uss;
    n= nss;
    m = mss;
    us = usss;
    f =fss;
    v = vss;
    theta = v/us;
    %theta = v/u;
    y = yss;
    c = css;
    lambda = lambdass;
    R = Rss;
    x = xss;
    J = Jss;
    Jprime = betta*Jss;
    w = wss;   
    d = x-w;
    REquity = ((1-delta)*J)/(J-d);
    sigma_z = sigma_zbar;
    RReal = R/Pii;
    ac = 0;    
    expGrowth = 0;
    S = exp(sss);
    SDF = betta;
    RFR = 1/SDF;
    Dividend = x*z - w;
 
    SDF1 = betta; SDF2 = betta; SDF3 = betta; SDF4 = betta; SDF5 = betta; SDF6 = betta; SDF7 = betta; SDF8 = betta;
    SDF9 = betta; SDF10 = betta; SDF11 = betta; SDF12 = betta;
    
    RE1 = 1/SDF; RE2 = 1/SDF; RE3 = 1/SDF; RE4 = 1/SDF; RE5 = 1/SDF; RE6 = 1/SDF; RE7 = 1/SDF; RE8 = 1/SDF;
    RE9 = 1/SDF; RE10 = 1/SDF; RE11 = 1/SDF; RE12 = 1/SDF;
    
    RFR1 = 1/SDF; RFR2 = 1/SDF; RFR3 = 1/SDF; RFR4 = 1/SDF; RFR5 = 1/SDF; RFR6 = 1/SDF; RFR7 = 1/SDF; RFR8 = 1/SDF;
    RFR9 = 1/SDF; RFR10 = 1/SDF; RFR11 = 1/SDF; RFR12 = 1/SDF;
    
    REcum = (1/SDF)^(12);
    RFRcum = (1/SDF)^(12);
    
    Ex1 = xss;   Ex2 =  xss;   Ex3 = xss;   Ex4 = xss;   Ex5 = xss;  Ex6 = xss;  Ex7 = xss;  Ex8 = xss;  Ex9 = xss;  Ex10 = xss; 
    Ec1 = css; Ec2 = css;Ec3 = css;Ec4 = css;Ec5 = css; Ec6 = css;Ec7 = css;Ec8 = css; Ec9 = css;Ec10 = css;
    Eu1 = uss;  Eu2 = uss;Eu3 = uss;Eu4 = uss;Eu5 = uss;Eu6 = uss;Eu7 = uss;Eu8 = uss;Eu9 = uss;Eu10 = uss;
    EJ1 = Jss; EJ2=Jss; EJ3=Jss; EJ4=Jss; EJ5=Jss; EJ6=Jss; EJ7=Jss; EJ8=Jss; EJ9=Jss; EJ10=Jss;
    Eac1=0; Eac2=0; Eac3=0; Eac4=0; Eac5=0; Eac6=0; Eac7=0; Eac8=0; Eac9=0; Eac10=0;
    Et1 = theta; Et2=theta; Et3=theta; Et4=theta; Et5=theta; Et6=theta; Et7=theta; Et8=theta; Et9=theta; Et10=theta;
    Ew1=w; Ew2=w; Ew3=w; Ew4=w; Ew5=w; Ew6=w; Ew7=w; Ew8=w; Ew9=w; Ew10=w;
    ERFR1=RFR; ERFR2=RFR; ERFR3=RFR; ERFR4=RFR; ERFR5=RFR; ERFR6=RFR; ERFR7=RFR; ERFR8=RFR; ERFR9=RFR; ERFR10=RFR;
    EPi1=Piiss; EPi2=Piiss; EPi3=Piiss; EPi4=Piiss; EPi5=Piiss; EPi6=Piiss; EPi7=Piiss; EPi8=Piiss; EPi9=Piiss; EPi10=Piiss;
    ER1=R; ER2=R; ER3=R; ER4=R; ER5=R; ER6=R; ER7=R; ER8=R; ER9=R;ER10=R;
    
    Elambda1 = lambda;
    CovJElambda1 = 0;
    JCheck = J;
    JNoCov = J;
    JNoElambda1 = J;
    
      
    sigma_R = sigma_RBar;
    RP = 0;
    RP1 = 0;
    
    n_Obs = n;
    z_Obs = z;
    sigma_z_Obs = sigma_z;
    logu = log(u);
    logtheta = log(theta);
    logz = log(z);

  
%% End own model equations

for iter = 1:length(M_.params) % Update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; % Auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);                % Get the steady state value of this variable
end
end