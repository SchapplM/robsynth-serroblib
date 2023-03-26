% Umwandlung der Kinematikparameter von S5RRRPR2V1 zu S5RRRPR2
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRRPR2V1
%   pkin_var=[a2 a3 d1 d2 d3 theta4]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRPR2
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d5 theta4]
%
% Siehe auch: S5RRRPR2_structural_kinematic_parameters.m
function pkin_gen = S5RRRPR2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 5, 6, 7, 9]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(8) = 0.0; % d5
