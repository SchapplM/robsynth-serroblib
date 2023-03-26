% Umwandlung der Kinematikparameter von S5RRRPR10V3 zu S5RRRPR10
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRRPR10V3
%   pkin_var=[a2 alpha2 d1 d2 theta4]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5RRRPR10
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d5 theta4]
%
% Siehe auch: S5RRRPR10_structural_kinematic_parameters.m
function pkin_gen = S5RRRPR10V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 5, 6, 7, 10]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(8) = 0.0; % d3
pkin_gen(9) = 0.0; % d5
