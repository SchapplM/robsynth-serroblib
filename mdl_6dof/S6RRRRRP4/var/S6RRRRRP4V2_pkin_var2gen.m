% Umwandlung der Kinematikparameter von S6RRRRRP4V2 zu S6RRRRRP4
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRRRP4V2
%   pkin_var=[a2 a3 a5 a6 d1 d2 d3 d5]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRRRP4
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d5]
%
% Siehe auch: S6RRRRRP4_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRP4V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 4, 5, 6, 7, 8, 10]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(9) = 0.0; % d4
