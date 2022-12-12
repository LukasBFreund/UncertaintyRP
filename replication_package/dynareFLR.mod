//=========================================================================
// Last updated: December 2022
//=========================================================================

//=========================================================================
//User Options
//=========================================================================
% Pre-processor variable: also compute GIRFs as per Andreasen (slower)
@# define OptionIRFAndreasen = 0

% Pre-processor variable to enable computation of expectations (slower)
@# define OptionExpectations = 0

% Pre-processor variable to enable computation of ergodic moments (slower)
@# define OptionMoments = 1

//--------------------------------------------------------------------------
%Pre-processor variable to switch between flexible (1) and sticky prices (0)
@# define OptionFlexprice = 0

% Pre-processor variable to use linearized NKPC (only usable with OptionFlexprice = 0)
@# define OptionPCLin = 1

% Pre-processor variable: habit formation on/off (1/0)
@# define OptionHabit = 1

% Pre-processor variable to switch to linearized marginal utility - equivalent to quad. utility
@# define OptionMULin = 0
   
% Pre-processor variable to turn off covariance term in JF
@# define OptionJNoCov = 0

% Pre-processor variable to linearize LoM for employment
@# define OptionLoMLin = 0

% Pre-processor variable to turn off E[lambda_{t+1}] effect
@# define OptionJNoElambda = 0

% Pre-processor variable to remove price adj cost from RC (default: 1)
@# define OptionNoPAC = 1

% Pre-processor variable to remove price adj cost from RC (default: 1)
@# define OptionNoVacCost = 1

//=========================================================================
//DECLARATION OF ENDOGENOUS VARIABLES
//=========================================================================
var
    c       ${c}$               (long_name='Consumption')
    R       ${R}$               (long_name='Short-term nominal interest rate (gross)')
    RReal   ${R^r}$             (long_name='Real Interest Rate')
    Pii     ${\Pi}$             (long_name='Inflation rate (gross)')
    y       ${y}$               (long_name='Aggregate output')
    x       ${x}$               (long_name='Relative price of intermediate goods')
    m       ${m}$               (long_name='Number of matches')                
    u       ${u}$               (long_name='Unemployment rate')
    us      ${u^s}$             (long_name='Number of searching workers')
    v       ${v}$               (long_name='Number of vacancies')
    f       ${f}$               (long_name='Job finding rate')
    h       ${h}$               (long_name='Vacancy filling rate')
    n       ${n}$               (long_name='Employment')
    J       ${J}$               (long_name='Value of a firm (with match)')
    w       ${w}$               (long_name='Actual wage')
    lambda  ${\lambda}$         (long_name='Marginal utility')
    theta   ${\theta}$          (long_name='Labor market tightness')   
    ac      ${ac}$              (long_name='Price adjustment costs')
    z       ${z}$               (long_name='Productivity')
    sigma_z ${\sigma_z}$        (long_name='Technology shock volatility')
    sigma_R ${\sigma_R}$        (long_name='Interest rate shock volatility')
    e_R     ${e_R}$             (long_name='Interest rate shock process')
    RP      ${RP}$              (long_name='Risk premium')
    S       ${S}$               (long_name='Consumption state')
    expGrowth ${expGrowth}$     (long_name='Expected growth')
    logu                        (long_name='Log of unemployment')
    logtheta                    (long_name='Log of tightness')
    logz                        (long_name='Log of productivity')
    RP1                         (long_name='Risk premium - 1 period ahead')
                   
// The following are auxililiary variables
    SDF     
    RFR             
    Dividend 
    Elambda1                    
    CovJElambda1                
    JCheck                      
    JNoCov
    JNoElambda1
    REquity 
    n_Obs 
    z_Obs 
    sigma_z_Obs


// Risk premium block
    SDF1 SDF2 SDF3 SDF4 SDF5 SDF6 SDF7 SDF8 SDF9 SDF10 SDF11 SDF12
    RE1 RE2 RE3 RE4 RE5 RE6 RE7 RE8 RE9 RE10 RE11 RE12
    RFR1 RFR2 RFR3 RFR4 RFR5 RFR6 RFR7 RFR8 RFR9 RFR10 RFR11 RFR12
    REcum RFRcum

// The following are expectations
    EJ1 
@#if OptionExpectations > 0
    Eu1 Eu2 Eu3 Eu4 Eu5 Eu6 Eu7 Eu8 Eu9 Eu10
    Ec1 Ec2 Ec3 Ec4 Ec5 Ec6 Ec7 Ec8 Ec9 Ec10 
        EJ2 EJ3 EJ4 EJ5 EJ6 EJ7 EJ8 EJ9 EJ10
    Eac1 Eac2 Eac3 Eac4 Eac5 Eac6 Eac7 Eac8 Eac9 Eac10
    Ex1 Ex2 Ex3 Ex4 Ex5 Ex6 Ex7 Ex8 Ex9 Ex10
    Et1 Et2 Et3 Et4 Et5 Et6 Et7 Et8 Et9 Et10
    Ew1 Ew2 Ew3 Ew4 Ew5 Ew6 Ew7 Ew8 Ew9 Ew10
    ERFR1 ERFR2 ERFR3 ERFR4 ERFR5 ERFR6 ERFR7 ERFR8 ERFR9 ERFR10
    EPi1 EPi2 EPi3 EPi4 EPi5 EPi6 EPi7 EPi8 EPi9 EPi10
    ER1 ER2 ER3 ER4 ER5 ER6 ER7 ER8 ER9 ER10
@#endif

;
 
//=========================================================================
//DECLARATION OF EXOGENOUS VARIABLES
//=========================================================================
varexo
    eps_z   ${\epsilon_z}$          (long_name='Technology innovation')
    eps_zu  ${\epsilon_{\sigma_z}}$ (long_name='Technology uncertainty innovation')
    eps_R   ${\epsilon_R}$          (long_name='Interest rate innovation')
    eps_Ru  ${\epsilon_{\sigma_R}}$ (long_name='Interest rate uncertainty innovation')
;
 
//=========================================================================
//DECLARATION OF PARAMETERS       
//=========================================================================
parameters
    betta       ${\betta}$              (long_name='Subjective discount factor')
    xi          ${\xi}$                 (long_name='Inverse elasticity of intertemporal substitution')
    eta         ${\eta}$                (long_name='Elasticity of substitution b/w differentiated products')
    alpha       ${\alpha}$              (long_name='Elasticity of matching w.r.t. u')
    delta       ${\rho}$                (long_name='Exogenous job destruction rate')      
    kappa       ${\kappa}$              (long_name='Vacancy posting fixed cost')
    Omega_p     ${\Omega_p}$            (long_name='Price adjustment cost')
    phi_pi      ${\phi_{\pi}}$          (long_name='Taylor rule parameter for inflation')
    phi_y       ${\phi_{y}}$            (long_name='Taylor rule parameter for output')
    phi_r       ${\phi_{r}}$            (long_name='Taylor rule interest rate smoothing')
    omega       ${\omega}$              (long_name='Bargaining weight for workers')
    rho_z       ${\rho_z}$              (long_name='Technology shock, persistence')
    rho_sigma_z ${\rho_{\sigma_z}}$     (long_name='Technoloy uncertainty shock, persistence')
    sigma_zbar  ${\bar{\sigma}_z}$      (long_name='Technology uncertainty shock, mean value')
    sigma_sigma_z ${\sigma_{\sigma_z}}$ (long_name='Technology uncertainty shock, sd of innovation')
    chi           ${\chi}$              (long_name='Strike value')
    psi           ${\bar{\mu}}$         (long_name='Matching efficiency')
    CCHabit1    ${\rho_{s}$             (long_name='Habit formation, persistence')
    CCHabit2    ${State}$               (long_name='Habit formation, average state')
    
// The following are additional shock process parameters
    sigma_RBar 
    rho_sigma_R
    sigma_sigma_R
    rho_e_R

// The following are calibrated steady-state values of variables
    zBar        ${\bar{Z}}$             (long_name='Technology shock, mean value')
    PiiBar      ${\bar{\pi}}$           (long_name='Inflation, steady-state value')
    uBar        ${\bar{u}}$             (long_name='Unemployment rate, steady-state value')                                                                  
    hBar        ${\bar{q}^v}$           (long_name='Vacancy filling rate, steady-state value')    

// The following are steady-state values
    css         ${\bar{C}}$
    Rss         ${\bar{R}}$
    Piiss       ${\bar{\pi}}$
    yss         ${\bar{Y}}$
    xss         ${\bar{x}}$
    mss         ${\bar{m}}$
    vss         ${\bar{v}}$
    fss         ${\bar{f}$
    hss         ${\bar{h}}$
    nss         ${\bar{n}}$
    uss         ${\bar{u}}$
    usss        ${\bar{u}^s}$
    Jss         ${\bar{J}^F}$
    lambdass    ${\bar{\lambda}}$
    zss         ${\bar{Z}}$               
    wss         ${\bar{w}}$
    sss         ${\bar{s}}$
    thetass     ${\bar{\theta}}$
    
// auxiliary (not used in default but can be commented in to get exact value for Omega_p instead of LL supplied value)
StickyDuration  ${\hat{\Omega}_p}$      (long_name='Duration of price stickiness')  
thetp           ${\tilde{\Omega}_p}$    (long_name='Calvo equivalent')

;

//=========================================================================
//PARAMETER VALUES                
//=========================================================================

load Parameters

% Parameters identical across specifications
set_param_value('betta',vPar.betta);
eta = 10;    
alpha = 0.5; 
PiiBar = 1.0016;                                    
hBar = 0.331;                         
uBar = 0.064; 
zBar = 1;  
delta = 0.028;                      
set_param_value('kappa',vPar.kappa);

% Habit
set_param_value('CCHabit1',vPar.CCHabit1);
set_param_value('CCHabit2',vPar.CCHabit2);

% Utility
set_param_value('xi',vPar.xi);

% Policy
set_param_value('phi_pi',vPar.phi_pi);
set_param_value('phi_y',vPar.phi_y);
set_param_value('phi_r',vPar.phi_r);

% Productivity shocks
set_param_value('rho_z',vPar.rho_z);
set_param_value('sigma_zbar',vPar.sigma_zbar);
set_param_value('sigma_sigma_z',vPar.sigma_sigma_z);
set_param_value('rho_sigma_z',vPar.rho_sigma_z);

% Interest rate shocks 
set_param_value('sigma_RBar',vPar.sigma_RBar);
set_param_value('rho_sigma_R',vPar.rho_sigma_R);
set_param_value('rho_e_R',vPar.rho_e_R);
set_param_value('sigma_sigma_R',vPar.sigma_sigma_R);


//--------------------------------------------------------------------------
//    SS Relationships
//--------------------------------------------------------------------------

zss = zBar;     
Piiss = PiiBar;
uss = uBar; 
hss = hBar;  
nss = 1-uss;
mss = delta*nss;
usss = 1-(1-delta)*nss;
sss = log(CCHabit2);
fss = mss/usss;
vss = mss/hss;
thetass= vss/usss;
yss = zss*nss;
psi = mss/(usss^alpha*vss^(1-alpha));

@#if OptionNoVacCost>0
css=yss;
@#else
css = yss-kappa*vss;
@#endif

@#if OptionMULin > 0
    @#if OptionHabit > 0
    lambdass = (CCHabit2)^(-xi)*(css - css*xi + css*xi)/css^(xi + 1);
    @#else
    lambdass = (CCHabit2)^(-xi)*(css - css*xi + css*xi)/css^(xi + 1);
    @#endif
@#else
    @#if OptionHabit > 0
    lambdass = (CCHabit2)^(-xi)*css^(-xi);
    @#else
    lambdass = (CCHabit2)^(-xi)*css^(-xi);
    @#endif
@#endif

Rss = Piiss/betta; 
xss = (eta-1)/eta;  

Jss = (1/hss)*(kappa);
wss = xss*zss-(1-betta*(1-delta))*Jss; 
set_param_value('chi',vPar.chi);
omega = (wss-chi)/(xss*zss-chi);

@#if OptionFlexprice > 0
StickyDuration = 0;        
thetp = 0;
Omega_p = 0;
@#else
StickyDuration = 9;         
thetp = 1-(1/StickyDuration); 
Omega_p = ((eta*xss/Piiss^2)*(thetp/(1-thetp)))/(1-thetp*betta);
@#endif

//=========================================================================
//MODEL EQUATIONS                 
//=========================================================================
model;
 
//--------------------------------------------------------------------------
//   Representative household
//--------------------------------------------------------------------------
[name = 'Habit']
expGrowth = logz(+1)-logz;
log(S) = (1-CCHabit1)*sss + CCHabit1*log(S(-1)) + ((1/CCHabit2)*(1-2*(log(S(-1))-sss))^(1/2)-1)*((logz-logz(-1))-expGrowth(-1));

[name = 'Marginal utility']
@#if OptionMULin > 0
    @#if OptionHabit > 0
    lambda = (S)^(-xi)*(css - c*xi + css*xi)/css^(xi + 1);
    @#else
    lambda = (CCHabit2)^(-xi)*(css - c*xi + css*xi)/css^(xi + 1);
    @#endif
@#else
    @#if OptionHabit > 0
    lambda = (S)^(-xi)*(c)^(-xi);
    @#else
    lambda = (CCHabit2)^(-xi)*(c)^(-xi);
    @#endif
@#endif

[name = 'Bond Euler']
1 = betta*R*((lambda(+1)/lambda)*(1/Pii(+1)));

//--------------------------------------------------------------------------
//   Labor market
//--------------------------------------------------------------------------

[name = 'Searching workers']
us = 1-(1-delta)*n(-1);

[name = 'Matching function']
m = psi*us^alpha*v^(1-alpha);

[name = 'Job finding rate ']
f= m/us;

[name = 'Vacancy filling rate']
h = m/v;

[name = 'Labor market tightness']
theta = v/us;


[name = 'Employment LoM']
@#if OptionLoMLin > 0
n = nss+(1-delta)*(n(-1)-nss)+fss*(us-usss)+usss*(f-fss); % linearized version (removes employment asymmetries given alpha = 0.5)
@#else
n = (1-delta)*n(-1)+f*us;   
@#endif
 
[name = 'Unemployment rate']
u = 1-n;

[name = 'Wage']
w = omega*(x*z)+(1-omega)*chi; 

//--------------------------------------------------------------------------
//  Intermediate goods sector
//--------------------------------------------------------------------------

[name = 'Production Function']
y = z*n;

[name = 'Free entry']
kappa = h*J;

[name = 'Equity Euler']
@#if OptionJNoCov > 0
@#if OptionJNoElambda > 0
J = x*z-w+(1-delta)*betta*(1/lambda)*(lambdass*J(+1)); % neither discounting nor covariance effects - but allowing for u'(c_t) to move
@#else
J = x*z-w+(1-delta)*betta*(1/lambda)*(Elambda1*J(+1)); % only discounting effects  - but allowing for u'(c_t) to move
@#endif
@#else
@#if OptionJNoElambda > 0
J = x*z-w+(1-delta)*betta*(1/lambda)*(lambdass*J(+1)+CovJElambda1); % only covariance effect  - but allowing for u'(c_t) to move
@#else
J = x*z-w+(1-delta)*betta*(lambda(+1)/lambda)*J(+1); % BASELINE
@#endif
@#endif

[name = 'Expected marginal utility']
Elambda1 = lambda(+1);
[name = 'Covariance term']
CovJElambda1 = (lambda(+1)-Elambda1)*(J(+1)-EJ1);
[name = 'Equity Euler - check']
JCheck = x*z-w+(1-delta)*betta*(1/lambda)*(Elambda1*EJ1+CovJElambda1);
[name = 'Equity Euler - no covariance']
JNoCov = x*z-w+(1-delta)*betta*(1/lambda)*(Elambda1*EJ1);
[name = 'Equity Euler - no expected marginal utility']
JNoElambda1 = x*z-w+(1-delta)*betta*(1/lambda)*(lambdass*EJ1+CovJElambda1);


//--------------------------------------------------------------------------
//   Retail sector
//--------------------------------------------------------------------------

[name = 'Phillips Curve']
@#if OptionPCLin > 0
Pii-Piiss=betta*(Pii(+1)-Piiss)+eta/(Omega_p*Piiss)*(x-xss);

@#else
x = ((eta-1)/eta) + (Omega_p/eta)*Pii*(Pii/Piiss-1) - betta*((lambda(+1)/lambda))*(y(+1)/y)*(Omega_p/eta)*Pii(+1)*(Pii(+1)/Piiss-1);
@#endif

//--------------------------------------------------------------------------
//  Aggregation
//--------------------------------------------------------------------------

[name = 'Price adjustment cost']
ac = (Omega_p/2)*(((Pii/Piiss-1))^2)*y;

[name = 'Resource constraint']
@#if OptionNoPAC > 0
@#if OptionNoVacCost > 0
y = c; 
@#else
y = c+kappa*v;
@#endif
@#else
@#if OptionNoVacCost > 0
y = c+ac; 
@#else
y = c+ac+kappa*v;
@#endif
@#endif

//--------------------------------------------------------------------------
//  Government Policy
//--------------------------------------------------------------------------
 
[name = 'Taylor rule']
log(R/STEADY_STATE(R))=phi_r*log(R(-1)/STEADY_STATE(R))+(1-phi_r)*(phi_pi*log(Pii/STEADY_STATE(Pii))+phi_y*log(y(+1)/STEADY_STATE(y)))+ e_R;

//--------------------------------------------------------------------------
//   Shock Processes
//--------------------------------------------------------------------------
 
[name = 'Technology shock']
z = (1-rho_z)* zBar+rho_z*z(-1)+sigma_z(-1)*eps_z;

[name = 'Technology uncertainty shock']
sigma_z=(1-rho_sigma_z)*sigma_zbar+rho_sigma_z*sigma_z(-1)+sigma_sigma_z*eps_zu;

[name = 'Interest rate shock']
e_R = rho_e_R*e_R(-1) + sigma_R(-1)*eps_R;

[name = 'Interest rate volatility']
sigma_R=(1-rho_sigma_R)*sigma_RBar+rho_sigma_R*sigma_R(-1)+sigma_sigma_R*eps_Ru;


//--------------------------------------------------------------------------
//   Additional Variables
//--------------------------------------------------------------------------

[name = 'Real interest rate']
RReal = R/Pii(+1);

[name = 'Return on equity']
REquity = ((1-delta)*J(+1))/(J-(x*z-w));

[name = 'Dividend']
Dividend = x*z - w;

[name = 'Stochastic discount factor']
SDF = betta*lambda(+1)/lambda;
SDF1 = betta*lambda(+1)/lambda;
SDF2 = betta*lambda(+2)/lambda(+1);
SDF3 = betta*lambda(+3)/lambda(+2);
SDF4 = betta*lambda(+4)/lambda(+3);
SDF5 = betta*lambda(+5)/lambda(+4);
SDF6 = betta*lambda(+6)/lambda(+5);
SDF7 = betta*lambda(+7)/lambda(+6);
SDF8 = betta*lambda(+8)/lambda(+7);
SDF9 = betta*lambda(+9)/lambda(+8);
SDF10 = betta*lambda(+10)/lambda(+9);
SDF11 = betta*lambda(+11)/lambda(+10);
SDF12 = betta*lambda(+12)/lambda(+11);


[name = 'Risk-free rate(s)']
RFR = 1/SDF;
RFR1 = 1/SDF1;
RFR2 = 1/SDF2;
RFR3 = 1/SDF3;
RFR4 = 1/SDF4;
RFR5 = 1/SDF5;
RFR6 = 1/SDF6;
RFR7 = 1/SDF7;
RFR8 = 1/SDF8;
RFR9 = 1/SDF9;
RFR10 = 1/SDF10;
RFR11 = 1/SDF11;
RFR12 = 1/SDF12;

RFRcum = RFR1*RFR2*RFR3*RFR4*RFR5*RFR6*RFR7*RFR8*RFR9*RFR10*RFR11*RFR12;

[name = 'Risk premium']
RP1 = RE1 - RFR;

RE1 = ((J(+1)*(1-delta))/(J-(x*z-w)));
RE2 = ((J(+2)*(1-delta))/(J(+1)-(x(+1)*z(+1)-w(+1))));
RE3 = ((J(+3)*(1-delta))/(J(+2)-(x(+2)*z(+2)-w(+2))));
RE4 = ((J(+4)*(1-delta))/(J(+3)-(x(+3)*z(+3)-w(+3))));
RE5 = ((J(+5)*(1-delta))/(J(+4)-(x(+4)*z(+4)-w(+4))));
RE6 = ((J(+6)*(1-delta))/(J(+5)-(x(+5)*z(+5)-w(+5))));
RE7 = ((J(+7)*(1-delta))/(J(+6)-(x(+6)*z(+6)-w(+6))));
RE8 = ((J(+8)*(1-delta))/(J(+7)-(x(+7)*z(+7)-w(+7))));
RE9 = ((J(+9)*(1-delta))/(J(+8)-(x(+8)*z(+8)-w(+8))));
RE10 = ((J(+10)*(1-delta))/(J(+9)-(x(+9)*z(+9)-w(+9))));
RE11 = ((J(+11)*(1-delta))/(J(+10)-(x(+10)*z(+10)-w(+10))));
RE12 = ((J(+12)*(1-delta))/(J(+11)-(x(+11)*z(+11)-w(+11))));

REcum = RE1*RE2*RE3*RE4*RE5*RE6*RE7*RE8*RE9*RE10*RE11*RE12;
RP = REcum - RFRcum; % 12-period ahead RP



//--------------------------------------------------------------------------
//   Expectations
//--------------------------------------------------------------------------

EJ1 = J(+1);

@#if OptionExpectations > 0

Ex1 = x(+1); Ex2 = x(+2); Ex3 = x(+3); Ex4 = x(+4); Ex5 = x(+5); Ex6=x(+6); Ex7=x(+7);Ex8=x(+8);Ex9=x(+9);Ex10=x(+10);
Ec1 = c(+1); Ec2 = c(+2); Ec3 = c(+3); Ec4 = c(+4); Ec5 = c(+5);Ec6=c(+6);Ec7=c(+7);Ec8=c(+8);Ec9=c(+9);Ec10=c(+10);
Eu1 = u(+1); Eu2 = u(+2); Eu3 = u(+3); Eu4 = u(+4); Eu5 = u(+5); Eu6=u(+6); Eu7=u(+7);Eu8=u(+8);Eu9=u(+9);Eu10=u(+10);
             EJ2 = J(+2); EJ3 = J(+3); EJ4 = J(+4); EJ5 = J(+5);EJ6=J(+6);EJ7=J(+7);EJ8=J(+8);EJ9=J(+9);EJ10=J(+10);
Eac1 = ac(+1); Eac2=ac(+2); Eac3=ac(+3); Eac4 = ac(+4); Eac5=ac(+5); Eac6=ac(+6); Eac7=ac(+7); Eac8=ac(+8); Eac9=ac(+9);Eac10=ac(+10);
Et1 = theta(+1); Et2=theta(+2); Et3=theta(+3); Et4=theta(+4); Et5=theta(+5); Et6=theta(+6); Et7=theta(+7); Et8=theta(+8); Et9=theta(+9); Et10=theta(+10);
Ew1=w(+1); Ew2=w(+2); Ew3=w(+3); Ew4=w(+4); Ew5=w(+5); Ew6=w(+6); Ew7=w(+7); Ew8=w(+8); Ew9=w(+9); Ew10=w(+10);
ERFR1=RFR(+1); ERFR2=RFR(+2); ERFR3=RFR(+3); ERFR4=RFR(+4); ERFR5=RFR(+5); ERFR6=RFR(+6); ERFR7=RFR(+7); ERFR8=RFR(+8); ERFR9=RFR(+9); ERFR10=RFR(+10);

EPi1=Pii(+1); EPi2=Pii(+2); EPi3=Pii(+3); EPi4=Pii(+4); EPi5=Pii(+5); EPi6=Pii(+6); EPi7=Pii(+7); EPi8=Pii(+8); EPi9=Pii(+9); EPi10=Pii(+10);
ER1=R(+1);ER2=R(+2); ER3=R(+3); ER4=R(+4); ER5=R(+5); ER6=R(+6); ER7=R(+7); ER8=R(+8); ER9=R(+9); ER10=R(+10);

@#endif


//--------------------------------------------------------------------------
//   Auxiliary
//--------------------------------------------------------------------------

// To trick Dynare into thinking that the following aren't state variables (for Andreasen method)
n_Obs = n;
z_Obs = z;
sigma_z_Obs = sigma_z;
logu = log(u);
logtheta = log(theta);
logz = log(z);

end;

//=========================================================================
//SHOCK VARIANCES                 
//=========================================================================
shocks;
var eps_z = 1;
var eps_zu = 1;           
var eps_R = 1;
var eps_Ru = 1;
end;

//=========================================================================
//STEADY STATE
//=========================================================================

% Provided in external steady-state file 

//--------------------------------------------------------------------------
// Check the starting values for the (det.) steady state
//--------------------------------------------------------------------------
resid;
 
//--------------------------------------------------------------------------
// Compute (det.) steady-state given the starting values
//--------------------------------------------------------------------------
steady;
 
//--------------------------------------------------------------------------
// Check Blanchard-Kahn-conditions
//--------------------------------------------------------------------------
check;


//=========================================================================
// POLICY FUNCTIONS
//=========================================================================

@#if OptionMoments > 0
stoch_simul(order=3,periods = 100000, drop = 50000,pruning,k_order_solver,irf=0) z logz c J u logu theta logtheta REquity RP RFR Dividend;
% hp_filter=100000
@#else
stoch_simul(order=3,pruning,k_order_solver,irf=0); 
@#endif

vNames = M_.endo_names; 

//=========================================================================
// Pure Uncertainty Shock IRF: Productivity
//=========================================================================
tic

% Specification
IrfPeriods = 60;                
BurninPeriods = 10000;   % periods for convergence from det. SS to EMAS
 
% Compute EMAS (Ergodic Mean in Absence of Shocks)
mEpsZero      = zeros(BurninPeriods+IrfPeriods,M_.exo_nbr); % shocks set to 0 to simulate without uncertainty
mIrfZero      = simult_(oo_.dr.ys,oo_.dr,mEpsZero,options_.order)'; % simulate series
vEMAS         = mIrfZero(1+BurninPeriods,:); % EMAS is any of the final points after burn-in
 
% Now simulate with an added impulse after the burn-in period
mEps = zeros(BurninPeriods+IrfPeriods,M_.exo_nbr);
mEps(1+BurninPeriods,strmatch('eps_zu',M_.exo_names,'exact')) = 1; % select productivity uncertainty

% Simulate model: starting from EMWS, burn-in to get us to EMAS, then impulse
mIrf = simult_(oo_.dr.ys,oo_.dr,mEps,options_.order)';
 
% Compute proportional deviation from EMAS
mIrfProp = (mIrf(1+BurninPeriods+1:1+BurninPeriods+IrfPeriods,:)-mIrfZero(1+BurninPeriods+1:1+BurninPeriods+IrfPeriods,:))./repmat(vEMAS,IrfPeriods,1); % only valid for variables not yet logged
 
% Adjust values that are virtually zero (not exactly due to use of a solver) to exactly zero to avoid confusing graphs
mIrfProp(abs(mIrfProp)<1e-12)=0;
 
mIRFProp_zUncertainty_EMAS = mIrfProp;

//=========================================================================
// Pure Uncertainty Shock IRF: Interest rate
//=========================================================================

% Specifications
IrfPeriods = 60;                
BurninPeriods = 10000;                 % periods for convergence from det. SS to EMAS
 
% Compute EMAS
mEpsZero      = zeros(BurninPeriods+IrfPeriods,M_.exo_nbr); % shocks set to 0 to simulate without uncertainty
mIrfZero      = simult_(oo_.dr.ys,oo_.dr,mEpsZero,options_.order)'; % simulate series
vEMAS         = mIrfZero(1+BurninPeriods,:); % EMAS is any of the final points after burn-in
 
% Now simulate with an added impulse after the burn-in period
mEps = zeros(BurninPeriods+IrfPeriods,M_.exo_nbr);
mEps(1+BurninPeriods,strmatch('eps_Ru',M_.exo_names,'exact')) = 1; % select productivity uncertainty

% Simulate model: starting from EMWS, burn-in to get us to EMAS, then impulse
mIrf = simult_(oo_.dr.ys,oo_.dr,mEps,options_.order)';
 
% Compute proportional deviation from EMAS
mIrfProp = (mIrf(1+BurninPeriods+1:1+BurninPeriods+IrfPeriods,:)-mIrfZero(1+BurninPeriods+1:1+BurninPeriods+IrfPeriods,:))./repmat(vEMAS,IrfPeriods,1); % only valid for variables not yet logged
 
% Adjust values that are virtually zero (not exactly due to use of a solver) to exactly zero to avoid confusing graphs
mIrfProp(abs(mIrfProp)<1e-12)=0;
 
mIRFProp_RUncertainty_EMAS = mIrfProp;


//=========================================================================
// GIRFs a la Andreasen et al. (2018)
//=========================================================================
@#if OptionIRFAndreasen > 0
tic

addpath('simAndMoments3order');

// Default step: a first-order approximation needed for RunDynarePruning
stoch_simul(order = 1, noprint, nomoments, irf = 0);
f_11 = [oo_.dr.ghx oo_.dr.ghu];
vNames = M_.endo_names;
vNames_cell = cellstr(vNames);
 
// Standard Dynare command 
stoch_simul(order = 3, noprint, nomoments, irf = 0);
 
// Options for running the pruning codes
optPruning.numSim         = 2000;                   
optPruning.seedNum        = 1; 
optPruning.orderApp       = options_.order;                                     
optPruning.computeIRF     = 1;
optPruning.ySelect     = vNames;
optPruning.plotIRF = 0;
optPruning.IRFlength = 60;
 
% Get a positive (+1) or negative (-1) shock
optPruning.shockSize = +1;
outDynare = RunDynarePruning(optPruning,oo_,M_,f_11);
vNames_Andreasen = outDynare.label_y;
 
vEMWS_Andreasen = outDynare.Mean_y';

% Recover IRFs and scale by EMAS
Order = 3;

mIRFProp_Andreasen_z_EMAS = outDynare.IRFy(:,:,Order,1)'./repmat(vEMAS,IrfPeriods,1);
mIRFProp_Andreasen_zVol_EMAS = outDynare.IRFy(:,:,Order,2)'./repmat(vEMAS,IrfPeriods,1);
mIRFProp_Andreasen_R_EMAS = outDynare.IRFy(:,:,Order,3)'./repmat(vEMAS,IrfPeriods,1);
mIRFProp_Andreasen_RVol_EMAS = outDynare.IRFy(:,:,Order,4)'./repmat(vEMAS,IrfPeriods,1);

% CAREFUL: the above only loads the IRFs for (what Dynare takes to be) controls! If you want to look at states
% you either have to declare them separately (as I've done e.g. for z, sigma_z, n) or use the following
mIRFProp_Andreasen_zVol_state = outDynare.IRFv(:,:,3,2)';

disp(['Run time for computation of IRFs a la Andreasen was ',num2str(toc),' seconds']);

@#endif

//=========================================================================
// SAVE STUFF
//=========================================================================

@#if OptionIRFAndreasen > 0
save("./Output/IRFs/IRFs_FLR_baseline", 'vEMAS', 'vNames', 'mIRFProp_zUncertainty_EMAS', 'mIRFProp_RUncertainty_EMAS', 'vNames_Andreasen', 'mIRFProp_Andreasen_R_EMAS', 'mIRFProp_Andreasen_RVol_EMAS', 'mIRFProp_Andreasen_zVol_EMAS')
@#else
save("./Output/IRFs/IRFs_FLR_baseline", 'vEMAS', 'vNames', 'mIRFProp_zUncertainty_EMAS', 'mIRFProp_RUncertainty_EMAS') 
@#endif 

