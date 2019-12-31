% Umwandlung der Kinematikparameter von S6RPRRRP8V1 zu S6RPRRRP8
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RPRRRP8V1
%   pkin_var=[a2 a3 a4 a6 d1 d3 d4]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPRRRP8
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 d5]
%
% Siehe auch: S6RPRRRP8_structural_kinematic_parameters.m
function pkin_gen = S6RPRRRP8V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 3, 5, 6, 7, 8]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d5
