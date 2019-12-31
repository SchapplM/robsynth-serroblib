% Umwandlung der Kinematikparameter von S6RRRRRR2V3 zu S6RRRRRR2
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRRRRR2V3
%   pkin_var=[a2 a3 a4 a6 d1 d2 d3 d4 d6]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRRRR2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d5 d6]
%
% Siehe auch: S6RRRRRR2_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRR2V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 3, 5, 6, 7, 8, 9, 11]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(10) = 0.0; % d5
