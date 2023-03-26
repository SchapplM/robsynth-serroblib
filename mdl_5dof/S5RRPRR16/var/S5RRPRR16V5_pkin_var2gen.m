% Umwandlung der Kinematikparameter von S5RRPRR16V5 zu S5RRPRR16
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRPRR16V5
%   pkin_var=[a2 a5 alpha2 d1 d2 d5]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPRR16
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d4 d5]
%
% Siehe auch: S5RRPRR16_structural_kinematic_parameters.m
function pkin_gen = S5RRPRR16V5_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 4, 5, 6, 7, 9]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(8) = 0.0; % d4
