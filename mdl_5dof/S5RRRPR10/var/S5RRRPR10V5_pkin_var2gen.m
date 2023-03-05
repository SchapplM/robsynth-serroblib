% Umwandlung der Kinematikparameter von S5RRRPR10V5 zu S5RRRPR10
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRRPR10V5
%   pkin_var=[a2 a3 alpha2 d1 d2 d3 theta4]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5RRRPR10
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d5 theta4]
%
% Siehe auch: S5RRRPR10_structural_kinematic_parameters.m
function pkin_gen = S5RRRPR10V5_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 5, 6, 7, 8, 10]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d5
