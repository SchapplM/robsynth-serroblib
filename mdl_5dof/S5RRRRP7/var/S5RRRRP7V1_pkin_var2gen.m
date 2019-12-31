% Umwandlung der Kinematikparameter von S5RRRRP7V1 zu S5RRRRP7
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRRRP7V1
%   pkin_var=[a3 a5 d1 d3]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRRP7
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4]
%
% Siehe auch: S5RRRRP7_structural_kinematic_parameters.m
function pkin_gen = S5RRRRP7V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([2, 4, 5, 7]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(3) = 0.0; % a4
pkin_gen(6) = 0.0; % d2
pkin_gen(8) = 0.0; % d4
