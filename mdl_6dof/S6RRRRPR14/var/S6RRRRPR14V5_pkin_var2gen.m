% Umwandlung der Kinematikparameter von S6RRRRPR14V5 zu S6RRRRPR14
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRPR14V5
%   pkin_var=[a4 a5 a6 d1 d4 d6 theta5]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRPR14
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d6 theta5]
%
% Siehe auch: S6RRRRPR14_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPR14V5_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([3, 4, 5, 8, 11, 12, 13]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(6) = pi/2; % alpha2
pkin_gen(7) = pi/2; % alpha3
pkin_gen(9) = 0.0; % d2
pkin_gen(10) = 0.0; % d3
