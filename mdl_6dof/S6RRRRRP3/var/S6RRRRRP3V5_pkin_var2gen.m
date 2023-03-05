% Umwandlung der Kinematikparameter von S6RRRRRP3V5 zu S6RRRRRP3
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRRRRP3V5
%   pkin_var=[a2 a3 a4 a5 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRRRP3
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d5]
%
% Siehe auch: S6RRRRRP3_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRP3V5_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 4, 6, 7, 8, 9, 10]) = pkin_var;

pkin_gen(5) = 0.0; % a6
