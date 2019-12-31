% Umwandlung der Kinematikparameter von S5RRPPP1V1 zu S5RRPPP1
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRPPP1V1
%   pkin_var=[a3 a4 a5 alpha3 d1 theta3]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRPPP1
%   pkin_gen=[a2 a3 a4 a5 alpha3 d1 d2 theta3]
%
% Siehe auch: S5RRPPP1_structural_kinematic_parameters.m
function pkin_gen = S5RRPPP1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([2, 3, 4, 5, 6, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(7) = 0.0; % d2
