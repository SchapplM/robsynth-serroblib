% Umwandlung der Kinematikparameter von S6RRPRRR13V5 zu S6RRPRRR13
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S6RRPRRR13V5
%   pkin_var=[a6 d1 d6]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR13
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 d6]
%
% Siehe auch: S6RRPRRR13_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRR13V5_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([5, 7, 11]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(6) = pi/2; % alpha2
pkin_gen(8) = 0.0; % d2
pkin_gen(9) = 0.0; % d4
pkin_gen(10) = 0.0; % d5
