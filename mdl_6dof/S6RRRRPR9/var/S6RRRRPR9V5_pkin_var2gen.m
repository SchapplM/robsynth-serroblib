% Umwandlung der Kinematikparameter von S6RRRRPR9V5 zu S6RRRRPR9
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S6RRRRPR9V5
%   pkin_var=[a4 d1 d4 theta5]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRPR9
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d4 d6 theta5]
%
% Siehe auch: S6RRRRPR9_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPR9V5_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([3, 7, 10, 12]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(6) = pi/2; % alpha2
pkin_gen(8) = 0.0; % d2
pkin_gen(9) = 0.0; % d3
pkin_gen(11) = 0.0; % d6
