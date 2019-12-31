% Umwandlung der Kinematikparameter von S5RRRPR4V1 zu S5RRRPR4
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRRPR4V1
%   pkin_var=[a2 a4 a5 d1 d2 d5]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRPR4
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d5]
%
% Siehe auch: S5RRRPR4_structural_kinematic_parameters.m
function pkin_gen = S5RRRPR4V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 3, 4, 5, 6, 8]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(7) = 0.0; % d3
