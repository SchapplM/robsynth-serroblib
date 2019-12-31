% Umwandlung der Kinematikparameter von S6PRRPPR1V2 zu S6PRRPPR1
% Eingabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6PRRPPR1V2
%   pkin_var=[a2 a4 a5 a6 alpha2 d2 d6 theta1 theta4 theta5]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRRPPR1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d3 d6 theta1 theta4 theta5]
%
% Siehe auch: S6PRRPPR1_structural_kinematic_parameters.m
function pkin_gen = S6PRRPPR1V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 3, 4, 5, 6, 7, 9, 10, 11, 12]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(8) = 0.0; % d3
