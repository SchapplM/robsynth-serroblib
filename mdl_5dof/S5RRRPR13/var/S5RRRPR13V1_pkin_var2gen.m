% Umwandlung der Kinematikparameter von S5RRRPR13V1 zu S5RRRPR13
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRRPR13V1
%   pkin_var=[a2 a4 a5 alpha2 d1 d2 d5]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRPR13
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d5]
%
% Siehe auch: S5RRRPR13_structural_kinematic_parameters.m
function pkin_gen = S5RRRPR13V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 3, 4, 5, 6, 7, 9]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(8) = 0.0; % d3
