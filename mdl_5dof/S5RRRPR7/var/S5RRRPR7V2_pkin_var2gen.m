% Umwandlung der Kinematikparameter von S5RRRPR7V2 zu S5RRRPR7
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRRPR7V2
%   pkin_var=[a3 d1 d3 theta4]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRPR7
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d5 theta4]
%
% Siehe auch: S5RRRPR7_structural_kinematic_parameters.m
function pkin_gen = S5RRRPR7V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([2, 5, 7, 9]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(6) = 0.0; % d2
pkin_gen(8) = 0.0; % d5
