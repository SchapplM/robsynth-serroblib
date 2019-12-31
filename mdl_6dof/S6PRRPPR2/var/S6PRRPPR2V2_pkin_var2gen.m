% Umwandlung der Kinematikparameter von S6PRRPPR2V2 zu S6PRRPPR2
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PRRPPR2V2
%   pkin_var=[a2 a4 a5 a6 alpha2 d2 d6 theta1 theta4]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6PRRPPR2
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d3 d6 theta1 theta4]
%
% Siehe auch: S6PRRPPR2_structural_kinematic_parameters.m
function pkin_gen = S6PRRPPR2V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 3, 4, 5, 6, 7, 9, 10, 11]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(8) = 0.0; % d3
