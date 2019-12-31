% Umwandlung der Kinematikparameter von S6RRRPPP1V1 zu S6RRRPPP1
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRPPP1V1
%   pkin_var=[a4 a5 a6 alpha4 d1 theta4]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRPPP1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha4 d1 d2 d3 theta4]
%
% Siehe auch: S6RRRPPP1_structural_kinematic_parameters.m
function pkin_gen = S6RRRPPP1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([3, 4, 5, 6, 7, 10]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(8) = 0.0; % d2
pkin_gen(9) = 0.0; % d3
