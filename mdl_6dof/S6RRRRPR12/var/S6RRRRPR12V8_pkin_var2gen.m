% Umwandlung der Kinematikparameter von S6RRRRPR12V8 zu S6RRRRPR12
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RRRRPR12V8
%   pkin_var=[a3 alpha3 d1 d3 theta5]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRPR12
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d6 theta5]
%
% Siehe auch: S6RRRRPR12_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPR12V8_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([2, 7, 8, 10, 13]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(6) = pi/2; % alpha2
pkin_gen(9) = 0.0; % d2
pkin_gen(11) = 0.0; % d4
pkin_gen(12) = 0.0; % d6
