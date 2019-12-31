% Umwandlung der Kinematikparameter von S6PPRPRR1V2 zu S6PPRPRR1
% Eingabe:
% pkin_var (11x1) double
%   Kinematikparameter (pkin) von S6PPRPRR1V2
%   pkin_var=[a2 a3 a4 a5 alpha2 alpha3 d3 d5 theta1 theta2 theta4]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6PPRPRR1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d3 d5 d6 theta1 theta2 theta4]
%
% Siehe auch: S6PPRPRR1_structural_kinematic_parameters.m
function pkin_gen = S6PPRPRR1V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13]) = pkin_var;

pkin_gen(5) = 0.0; % a6
pkin_gen(10) = 0.0; % d6
