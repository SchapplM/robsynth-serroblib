% Umwandlung der Kinematikparameter von S6RRRRRP2V2 zu S6RRRRRP2
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRRRP2V2
%   pkin_var=[a2 a3 a4 a6 d1 d2 d3 d4]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRRRP2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d5]
%
% Siehe auch: S6RRRRRP2_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRP2V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 5, 6, 7, 8, 9]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(10) = 0.0; % d5
