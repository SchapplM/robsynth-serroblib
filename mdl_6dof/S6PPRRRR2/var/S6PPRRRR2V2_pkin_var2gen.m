% Umwandlung der Kinematikparameter von S6PPRRRR2V2 zu S6PPRRRR2
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PPRRRR2V2
%   pkin_var=[a2 a3 a6 alpha2 alpha3 d3 d6 theta1 theta2]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6PPRRRR2
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d3 d4 d5 d6 theta1 theta2]
%
% Siehe auch: S6PPRRRR2_structural_kinematic_parameters.m
function pkin_gen = S6PPRRRR2V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([1, 2, 5, 6, 7, 8, 11, 12, 13]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d4
pkin_gen(10) = 0.0; % d5
