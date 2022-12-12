function [residual, g1, g2] = dynareFLR_noRP_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 102, 1);

%
% Model equations
%

T3 = (-1);
T27 = 1/params(19)*(1-2*(log(y(24))-params(44)))^0.5-1;
T100 = params(1)*(1-params(5))*1/y(16);
T147 = params(7)/2*(y(4)/params(30)-1)^2;
T219 = (1-params(5))*y(14)/(y(14)-(y(6)*y(19)-y(15)));
T413 = ((y(2))-y(2))/((y(2))*(y(2)))/(y(2)/(y(2)));
T454 = (-((-((1-params(5))*y(14)*(-y(19))))/((y(14)-(y(6)*y(19)-y(15)))*(y(14)-(y(6)*y(19)-y(15))))));
T490 = (-(((1-params(5))*(y(14)-(y(6)*y(19)-y(15)))-(1-params(5))*y(14))/((y(14)-(y(6)*y(19)-y(15)))*(y(14)-(y(6)*y(19)-y(15))))));
T493 = (-((-((1-params(5))*y(14)))/((y(14)-(y(6)*y(19)-y(15)))*(y(14)-(y(6)*y(19)-y(15))))));
T496 = params(1)*(1-params(5))*T3/(y(16)*y(16));
T515 = (-((-((1-params(5))*y(14)*(-y(6))))/((y(14)-(y(6)*y(19)-y(15)))*(y(14)-(y(6)*y(19)-y(15))))));
residual(1) = y(25);
lhs =log(y(24));
rhs =(1-params(18))*params(44)+log(y(24))*params(18)+T27*(-y(25));
residual(2)= lhs-rhs;
lhs =y(16);
rhs =y(24)^(-params(2))*y(1)^(-params(2));
residual(3)= lhs-rhs;
lhs =1;
rhs =params(1)*y(2)*1/y(4);
residual(4)= lhs-rhs;
lhs =y(9);
rhs =1-(1-params(5))*y(13);
residual(5)= lhs-rhs;
lhs =y(7);
rhs =params(17)*y(9)^params(4)*y(10)^(1-params(4));
residual(6)= lhs-rhs;
lhs =y(11);
rhs =y(7)/y(9);
residual(7)= lhs-rhs;
lhs =y(12);
rhs =y(7)/y(10);
residual(8)= lhs-rhs;
lhs =y(17);
rhs =y(10)/y(9);
residual(9)= lhs-rhs;
lhs =y(13);
rhs =(1-params(5))*y(13)+y(9)*y(11);
residual(10)= lhs-rhs;
lhs =y(8);
rhs =1-y(13);
residual(11)= lhs-rhs;
lhs =y(15);
rhs =params(11)*y(6)*y(19)+(1-params(11))*params(16);
residual(12)= lhs-rhs;
lhs =y(5);
rhs =y(13)*y(19);
residual(13)= lhs-rhs;
lhs =params(6);
rhs =y(12)*y(14);
residual(14)= lhs-rhs;
lhs =y(14);
rhs =y(6)*y(19)-y(15)+T100*y(14)*y(33);
residual(15)= lhs-rhs;
lhs =y(33);
rhs =y(16);
residual(16)= lhs-rhs;
lhs =y(34);
rhs =(y(16)-y(33))*(y(14)-y(80));
residual(17)= lhs-rhs;
lhs =y(35);
rhs =y(6)*y(19)-y(15)+T100*(y(34)+y(33)*y(80));
residual(18)= lhs-rhs;
lhs =y(36);
rhs =y(6)*y(19)-y(15)+T100*y(33)*y(80);
residual(19)= lhs-rhs;
lhs =y(37);
rhs =y(6)*y(19)-y(15)+T100*(y(34)+y(80)*params(41));
residual(20)= lhs-rhs;
lhs =y(4)-params(30);
rhs =params(1)*(y(4)-params(30))+params(3)/(params(30)*params(7))*(y(6)-params(32));
residual(21)= lhs-rhs;
lhs =y(18);
rhs =y(5)*T147;
residual(22)= lhs-rhs;
lhs =y(5);
rhs =y(1);
residual(23)= lhs-rhs;
lhs =log(y(2)/(y(2)));
rhs =log(y(2)/(y(2)))*params(10)+(1-params(10))*(params(8)*log(y(4)/(y(4)))+params(9)*log(y(5)/(y(5))))+y(22);
residual(24)= lhs-rhs;
lhs =y(19);
rhs =(1-params(12))*params(24)+y(19)*params(12)+y(20)*x(1);
residual(25)= lhs-rhs;
lhs =y(20);
rhs =(1-params(13))*params(14)+y(20)*params(13)+params(15)*x(2);
residual(26)= lhs-rhs;
lhs =y(22);
rhs =y(22)*params(23)+y(21)*x(3);
residual(27)= lhs-rhs;
lhs =y(21);
rhs =(1-params(21))*params(20)+y(21)*params(21)+params(22)*x(4);
residual(28)= lhs-rhs;
lhs =y(3);
rhs =y(2)/y(4);
residual(29)= lhs-rhs;
lhs =y(38);
rhs =T219;
residual(30)= lhs-rhs;
lhs =y(32);
rhs =y(6)*y(19)-y(15);
residual(31)= lhs-rhs;
lhs =y(30);
rhs =y(16)*params(1)/y(16);
residual(32)= lhs-rhs;
lhs =y(42);
rhs =y(16)*params(1)/y(16);
residual(33)= lhs-rhs;
lhs =y(43);
rhs =y(81);
residual(34)= lhs-rhs;
lhs =y(44);
rhs =y(82);
residual(35)= lhs-rhs;
lhs =y(45);
rhs =y(83);
residual(36)= lhs-rhs;
lhs =y(46);
rhs =y(84);
residual(37)= lhs-rhs;
lhs =y(47);
rhs =y(85);
residual(38)= lhs-rhs;
lhs =y(48);
rhs =y(86);
residual(39)= lhs-rhs;
lhs =y(49);
rhs =y(87);
residual(40)= lhs-rhs;
lhs =y(50);
rhs =y(88);
residual(41)= lhs-rhs;
lhs =y(51);
rhs =y(89);
residual(42)= lhs-rhs;
lhs =y(52);
rhs =y(90);
residual(43)= lhs-rhs;
lhs =y(53);
rhs =y(91);
residual(44)= lhs-rhs;
lhs =y(31);
rhs =1/y(30);
residual(45)= lhs-rhs;
lhs =y(66);
rhs =1/y(42);
residual(46)= lhs-rhs;
lhs =y(67);
rhs =1/y(43);
residual(47)= lhs-rhs;
lhs =y(68);
rhs =1/y(44);
residual(48)= lhs-rhs;
lhs =y(69);
rhs =1/y(45);
residual(49)= lhs-rhs;
lhs =y(70);
rhs =1/y(46);
residual(50)= lhs-rhs;
lhs =y(71);
rhs =1/y(47);
residual(51)= lhs-rhs;
lhs =y(72);
rhs =1/y(48);
residual(52)= lhs-rhs;
lhs =y(73);
rhs =1/y(49);
residual(53)= lhs-rhs;
lhs =y(74);
rhs =1/y(50);
residual(54)= lhs-rhs;
lhs =y(75);
rhs =1/y(51);
residual(55)= lhs-rhs;
lhs =y(76);
rhs =1/y(52);
residual(56)= lhs-rhs;
lhs =y(77);
rhs =1/y(53);
residual(57)= lhs-rhs;
lhs =y(79);
rhs =y(77)*y(76)*y(75)*y(74)*y(73)*y(72)*y(71)*y(70)*y(69)*y(68)*y(66)*y(67);
residual(58)= lhs-rhs;
lhs =y(29);
rhs =y(54)-y(31);
residual(59)= lhs-rhs;
lhs =y(54);
rhs =T219;
residual(60)= lhs-rhs;
lhs =y(55);
rhs =y(92);
residual(61)= lhs-rhs;
lhs =y(56);
rhs =y(93);
residual(62)= lhs-rhs;
lhs =y(57);
rhs =y(94);
residual(63)= lhs-rhs;
lhs =y(58);
rhs =y(95);
residual(64)= lhs-rhs;
lhs =y(59);
rhs =y(96);
residual(65)= lhs-rhs;
lhs =y(60);
rhs =y(97);
residual(66)= lhs-rhs;
lhs =y(61);
rhs =y(98);
residual(67)= lhs-rhs;
lhs =y(62);
rhs =y(99);
residual(68)= lhs-rhs;
lhs =y(63);
rhs =y(100);
residual(69)= lhs-rhs;
lhs =y(64);
rhs =y(101);
residual(70)= lhs-rhs;
lhs =y(65);
rhs =y(102);
residual(71)= lhs-rhs;
lhs =y(78);
rhs =y(65)*y(64)*y(63)*y(62)*y(61)*y(60)*y(59)*y(58)*y(57)*y(56)*y(54)*y(55);
residual(72)= lhs-rhs;
lhs =y(23);
rhs =y(78)-y(79);
residual(73)= lhs-rhs;
lhs =y(80);
rhs =y(14);
residual(74)= lhs-rhs;
lhs =y(39);
rhs =y(13);
residual(75)= lhs-rhs;
lhs =y(40);
rhs =y(19);
residual(76)= lhs-rhs;
lhs =y(41);
rhs =y(20);
residual(77)= lhs-rhs;
lhs =y(26);
rhs =log(y(8));
residual(78)= lhs-rhs;
lhs =y(27);
rhs =log(y(17));
residual(79)= lhs-rhs;
lhs =y(28);
rhs =log(y(19));
residual(80)= lhs-rhs;
lhs =y(81);
rhs =y(16)*params(1)/y(16);
residual(81)= lhs-rhs;
lhs =y(82);
rhs =y(81);
residual(82)= lhs-rhs;
lhs =y(83);
rhs =y(82);
residual(83)= lhs-rhs;
lhs =y(84);
rhs =y(83);
residual(84)= lhs-rhs;
lhs =y(85);
rhs =y(84);
residual(85)= lhs-rhs;
lhs =y(86);
rhs =y(85);
residual(86)= lhs-rhs;
lhs =y(87);
rhs =y(86);
residual(87)= lhs-rhs;
lhs =y(88);
rhs =y(87);
residual(88)= lhs-rhs;
lhs =y(89);
rhs =y(88);
residual(89)= lhs-rhs;
lhs =y(90);
rhs =y(89);
residual(90)= lhs-rhs;
lhs =y(91);
rhs =y(90);
residual(91)= lhs-rhs;
lhs =y(92);
rhs =T219;
residual(92)= lhs-rhs;
lhs =y(93);
rhs =y(92);
residual(93)= lhs-rhs;
lhs =y(94);
rhs =y(93);
residual(94)= lhs-rhs;
lhs =y(95);
rhs =y(94);
residual(95)= lhs-rhs;
lhs =y(96);
rhs =y(95);
residual(96)= lhs-rhs;
lhs =y(97);
rhs =y(96);
residual(97)= lhs-rhs;
lhs =y(98);
rhs =y(97);
residual(98)= lhs-rhs;
lhs =y(99);
rhs =y(98);
residual(99)= lhs-rhs;
lhs =y(100);
rhs =y(99);
residual(100)= lhs-rhs;
lhs =y(101);
rhs =y(100);
residual(101)= lhs-rhs;
lhs =y(102);
rhs =y(101);
residual(102)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(102, 102);

  %
  % Jacobian matrix
  %

  g1(1,25)=1;
  g1(2,24)=1/y(24)-(params(18)*1/y(24)+(-y(25))*1/params(19)*(-(2*1/y(24)))*0.5*(1-2*(log(y(24))-params(44)))^(-0.5));
  g1(2,25)=T27;
  g1(3,1)=(-(y(24)^(-params(2))*getPowerDeriv(y(1),(-params(2)),1)));
  g1(3,16)=1;
  g1(3,24)=(-(y(1)^(-params(2))*getPowerDeriv(y(24),(-params(2)),1)));
  g1(4,2)=(-(params(1)*1/y(4)));
  g1(4,4)=(-(params(1)*y(2)*T3/(y(4)*y(4))));
  g1(5,9)=1;
  g1(5,13)=1-params(5);
  g1(6,7)=1;
  g1(6,9)=(-(y(10)^(1-params(4))*params(17)*getPowerDeriv(y(9),params(4),1)));
  g1(6,10)=(-(params(17)*y(9)^params(4)*getPowerDeriv(y(10),1-params(4),1)));
  g1(7,7)=(-(1/y(9)));
  g1(7,9)=(-((-y(7))/(y(9)*y(9))));
  g1(7,11)=1;
  g1(8,7)=(-(1/y(10)));
  g1(8,10)=(-((-y(7))/(y(10)*y(10))));
  g1(8,12)=1;
  g1(9,9)=(-((-y(10))/(y(9)*y(9))));
  g1(9,10)=(-(1/y(9)));
  g1(9,17)=1;
  g1(10,9)=(-y(11));
  g1(10,11)=(-y(9));
  g1(10,13)=1-(1-params(5));
  g1(11,8)=1;
  g1(11,13)=1;
  g1(12,6)=(-(params(11)*y(19)));
  g1(12,15)=1;
  g1(12,19)=(-(params(11)*y(6)));
  g1(13,5)=1;
  g1(13,13)=(-y(19));
  g1(13,19)=(-y(13));
  g1(14,12)=(-y(14));
  g1(14,14)=(-y(12));
  g1(15,6)=(-y(19));
  g1(15,14)=1-T100*y(33);
  g1(15,15)=1;
  g1(15,16)=(-(y(14)*y(33)*T496));
  g1(15,19)=(-y(6));
  g1(15,33)=(-(y(14)*T100));
  g1(16,16)=T3;
  g1(16,33)=1;
  g1(17,14)=(-(y(16)-y(33)));
  g1(17,16)=(-(y(14)-y(80)));
  g1(17,33)=y(14)-y(80);
  g1(17,34)=1;
  g1(17,80)=y(16)-y(33);
  g1(18,6)=(-y(19));
  g1(18,15)=1;
  g1(18,16)=(-((y(34)+y(33)*y(80))*T496));
  g1(18,19)=(-y(6));
  g1(18,33)=(-(T100*y(80)));
  g1(18,34)=(-T100);
  g1(18,35)=1;
  g1(18,80)=(-(T100*y(33)));
  g1(19,6)=(-y(19));
  g1(19,15)=1;
  g1(19,16)=(-(y(33)*y(80)*T496));
  g1(19,19)=(-y(6));
  g1(19,33)=(-(T100*y(80)));
  g1(19,36)=1;
  g1(19,80)=(-(T100*y(33)));
  g1(20,6)=(-y(19));
  g1(20,15)=1;
  g1(20,16)=(-((y(34)+y(80)*params(41))*T496));
  g1(20,19)=(-y(6));
  g1(20,34)=(-T100);
  g1(20,37)=1;
  g1(20,80)=(-(T100*params(41)));
  g1(21,4)=1-params(1);
  g1(21,6)=(-(params(3)/(params(30)*params(7))));
  g1(22,4)=(-(y(5)*params(7)/2*1/params(30)*2*(y(4)/params(30)-1)));
  g1(22,5)=(-T147);
  g1(22,18)=1;
  g1(23,1)=T3;
  g1(23,5)=1;
  g1(24,2)=T413-params(10)*T413;
  g1(24,4)=(-((1-params(10))*params(8)*((y(4))-y(4))/((y(4))*(y(4)))/(y(4)/(y(4)))));
  g1(24,5)=(-((1-params(10))*params(9)*((y(5))-y(5))/((y(5))*(y(5)))/(y(5)/(y(5)))));
  g1(24,22)=T3;
  g1(25,19)=1-params(12);
  g1(25,20)=(-x(1));
  g1(26,20)=1-params(13);
  g1(27,21)=(-x(3));
  g1(27,22)=1-params(23);
  g1(28,21)=1-params(21);
  g1(29,2)=(-(1/y(4)));
  g1(29,3)=1;
  g1(29,4)=(-((-y(2))/(y(4)*y(4))));
  g1(30,6)=T454;
  g1(30,14)=T490;
  g1(30,15)=T493;
  g1(30,19)=T515;
  g1(30,38)=1;
  g1(31,6)=(-y(19));
  g1(31,15)=1;
  g1(31,19)=(-y(6));
  g1(31,32)=1;
  g1(32,30)=1;
  g1(33,42)=1;
  g1(34,43)=1;
  g1(34,81)=T3;
  g1(35,44)=1;
  g1(35,82)=T3;
  g1(36,45)=1;
  g1(36,83)=T3;
  g1(37,46)=1;
  g1(37,84)=T3;
  g1(38,47)=1;
  g1(38,85)=T3;
  g1(39,48)=1;
  g1(39,86)=T3;
  g1(40,49)=1;
  g1(40,87)=T3;
  g1(41,50)=1;
  g1(41,88)=T3;
  g1(42,51)=1;
  g1(42,89)=T3;
  g1(43,52)=1;
  g1(43,90)=T3;
  g1(44,53)=1;
  g1(44,91)=T3;
  g1(45,30)=(-(T3/(y(30)*y(30))));
  g1(45,31)=1;
  g1(46,42)=(-(T3/(y(42)*y(42))));
  g1(46,66)=1;
  g1(47,43)=(-(T3/(y(43)*y(43))));
  g1(47,67)=1;
  g1(48,44)=(-(T3/(y(44)*y(44))));
  g1(48,68)=1;
  g1(49,45)=(-(T3/(y(45)*y(45))));
  g1(49,69)=1;
  g1(50,46)=(-(T3/(y(46)*y(46))));
  g1(50,70)=1;
  g1(51,47)=(-(T3/(y(47)*y(47))));
  g1(51,71)=1;
  g1(52,48)=(-(T3/(y(48)*y(48))));
  g1(52,72)=1;
  g1(53,49)=(-(T3/(y(49)*y(49))));
  g1(53,73)=1;
  g1(54,50)=(-(T3/(y(50)*y(50))));
  g1(54,74)=1;
  g1(55,51)=(-(T3/(y(51)*y(51))));
  g1(55,75)=1;
  g1(56,52)=(-(T3/(y(52)*y(52))));
  g1(56,76)=1;
  g1(57,53)=(-(T3/(y(53)*y(53))));
  g1(57,77)=1;
  g1(58,66)=(-(y(77)*y(76)*y(75)*y(74)*y(73)*y(72)*y(71)*y(70)*y(69)*y(67)*y(68)));
  g1(58,67)=(-(y(77)*y(76)*y(75)*y(74)*y(73)*y(72)*y(71)*y(70)*y(69)*y(66)*y(68)));
  g1(58,68)=(-(y(77)*y(76)*y(75)*y(74)*y(73)*y(72)*y(71)*y(70)*y(69)*y(66)*y(67)));
  g1(58,69)=(-(y(77)*y(76)*y(75)*y(74)*y(73)*y(72)*y(71)*y(70)*y(68)*y(66)*y(67)));
  g1(58,70)=(-(y(77)*y(76)*y(75)*y(74)*y(73)*y(72)*y(71)*y(69)*y(68)*y(66)*y(67)));
  g1(58,71)=(-(y(77)*y(76)*y(75)*y(74)*y(73)*y(72)*y(70)*y(69)*y(68)*y(66)*y(67)));
  g1(58,72)=(-(y(77)*y(76)*y(75)*y(74)*y(73)*y(71)*y(70)*y(69)*y(68)*y(66)*y(67)));
  g1(58,73)=(-(y(77)*y(76)*y(75)*y(74)*y(72)*y(71)*y(70)*y(69)*y(68)*y(66)*y(67)));
  g1(58,74)=(-(y(77)*y(76)*y(75)*y(73)*y(72)*y(71)*y(70)*y(69)*y(68)*y(66)*y(67)));
  g1(58,75)=(-(y(77)*y(76)*y(74)*y(73)*y(72)*y(71)*y(70)*y(69)*y(68)*y(66)*y(67)));
  g1(58,76)=(-(y(77)*y(75)*y(74)*y(73)*y(72)*y(71)*y(70)*y(69)*y(68)*y(66)*y(67)));
  g1(58,77)=(-(y(76)*y(75)*y(74)*y(73)*y(72)*y(71)*y(70)*y(69)*y(68)*y(66)*y(67)));
  g1(58,79)=1;
  g1(59,29)=1;
  g1(59,31)=1;
  g1(59,54)=T3;
  g1(60,6)=T454;
  g1(60,14)=T490;
  g1(60,15)=T493;
  g1(60,19)=T515;
  g1(60,54)=1;
  g1(61,55)=1;
  g1(61,92)=T3;
  g1(62,56)=1;
  g1(62,93)=T3;
  g1(63,57)=1;
  g1(63,94)=T3;
  g1(64,58)=1;
  g1(64,95)=T3;
  g1(65,59)=1;
  g1(65,96)=T3;
  g1(66,60)=1;
  g1(66,97)=T3;
  g1(67,61)=1;
  g1(67,98)=T3;
  g1(68,62)=1;
  g1(68,99)=T3;
  g1(69,63)=1;
  g1(69,100)=T3;
  g1(70,64)=1;
  g1(70,101)=T3;
  g1(71,65)=1;
  g1(71,102)=T3;
  g1(72,54)=(-(y(65)*y(64)*y(63)*y(62)*y(61)*y(60)*y(59)*y(58)*y(57)*y(55)*y(56)));
  g1(72,55)=(-(y(65)*y(64)*y(63)*y(62)*y(61)*y(60)*y(59)*y(58)*y(57)*y(54)*y(56)));
  g1(72,56)=(-(y(65)*y(64)*y(63)*y(62)*y(61)*y(60)*y(59)*y(58)*y(57)*y(54)*y(55)));
  g1(72,57)=(-(y(65)*y(64)*y(63)*y(62)*y(61)*y(60)*y(59)*y(58)*y(56)*y(54)*y(55)));
  g1(72,58)=(-(y(65)*y(64)*y(63)*y(62)*y(61)*y(60)*y(59)*y(57)*y(56)*y(54)*y(55)));
  g1(72,59)=(-(y(65)*y(64)*y(63)*y(62)*y(61)*y(60)*y(58)*y(57)*y(56)*y(54)*y(55)));
  g1(72,60)=(-(y(65)*y(64)*y(63)*y(62)*y(61)*y(59)*y(58)*y(57)*y(56)*y(54)*y(55)));
  g1(72,61)=(-(y(65)*y(64)*y(63)*y(62)*y(60)*y(59)*y(58)*y(57)*y(56)*y(54)*y(55)));
  g1(72,62)=(-(y(65)*y(64)*y(63)*y(61)*y(60)*y(59)*y(58)*y(57)*y(56)*y(54)*y(55)));
  g1(72,63)=(-(y(65)*y(64)*y(62)*y(61)*y(60)*y(59)*y(58)*y(57)*y(56)*y(54)*y(55)));
  g1(72,64)=(-(y(65)*y(63)*y(62)*y(61)*y(60)*y(59)*y(58)*y(57)*y(56)*y(54)*y(55)));
  g1(72,65)=(-(y(64)*y(63)*y(62)*y(61)*y(60)*y(59)*y(58)*y(57)*y(56)*y(54)*y(55)));
  g1(72,78)=1;
  g1(73,23)=1;
  g1(73,78)=T3;
  g1(73,79)=1;
  g1(74,14)=T3;
  g1(74,80)=1;
  g1(75,13)=T3;
  g1(75,39)=1;
  g1(76,19)=T3;
  g1(76,40)=1;
  g1(77,20)=T3;
  g1(77,41)=1;
  g1(78,8)=(-(1/y(8)));
  g1(78,26)=1;
  g1(79,17)=(-(1/y(17)));
  g1(79,27)=1;
  g1(80,19)=(-(1/y(19)));
  g1(80,28)=1;
  g1(81,81)=1;
  g1(82,81)=T3;
  g1(82,82)=1;
  g1(83,82)=T3;
  g1(83,83)=1;
  g1(84,83)=T3;
  g1(84,84)=1;
  g1(85,84)=T3;
  g1(85,85)=1;
  g1(86,85)=T3;
  g1(86,86)=1;
  g1(87,86)=T3;
  g1(87,87)=1;
  g1(88,87)=T3;
  g1(88,88)=1;
  g1(89,88)=T3;
  g1(89,89)=1;
  g1(90,89)=T3;
  g1(90,90)=1;
  g1(91,90)=T3;
  g1(91,91)=1;
  g1(92,6)=T454;
  g1(92,14)=T490;
  g1(92,15)=T493;
  g1(92,19)=T515;
  g1(92,92)=1;
  g1(93,92)=T3;
  g1(93,93)=1;
  g1(94,93)=T3;
  g1(94,94)=1;
  g1(95,94)=T3;
  g1(95,95)=1;
  g1(96,95)=T3;
  g1(96,96)=1;
  g1(97,96)=T3;
  g1(97,97)=1;
  g1(98,97)=T3;
  g1(98,98)=1;
  g1(99,98)=T3;
  g1(99,99)=1;
  g1(100,99)=T3;
  g1(100,100)=1;
  g1(101,100)=T3;
  g1(101,101)=1;
  g1(102,101)=T3;
  g1(102,102)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],102,10404);
end
end