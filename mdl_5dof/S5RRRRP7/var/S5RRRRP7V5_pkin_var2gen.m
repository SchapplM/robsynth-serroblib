% Umwandlung der Kinematikparameter von S5RRRRP7V5 zu S5RRRRP7
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRRRP7V5
%   pkin_var=[a2 a3 a4 d1 d2 d3 d4]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRRP7
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4]
%
% Siehe auch: S5RRRRP7_structural_kinematic_parameters.m
function pkin_gen = S5RRRRP7V5_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 2, 3, 5, 6, 7, 8]) = pkin_var;

pkin_gen(4) = 0.0; % a5
