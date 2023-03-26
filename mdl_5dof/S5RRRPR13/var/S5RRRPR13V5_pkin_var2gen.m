% Umwandlung der Kinematikparameter von S5RRRPR13V5 zu S5RRRPR13
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRRPR13V5
%   pkin_var=[a2 a3 alpha2 d1 d2 d3]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRPR13
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d5]
%
% Siehe auch: S5RRRPR13_structural_kinematic_parameters.m
function pkin_gen = S5RRRPR13V5_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 5, 6, 7, 8]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d5
