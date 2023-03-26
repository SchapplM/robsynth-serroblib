% Umwandlung der Kinematikparameter von S6RRRRPR15V6 zu S6RRRRPR15
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRPR15V6
%   pkin_var=[a2 a3 alpha2 alpha3 d1 d2 d3]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRPR15
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d6]
%
% Siehe auch: S6RRRRPR15_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPR15V6_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 6, 7, 8, 9, 10]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(11) = 0.0; % d4
pkin_gen(12) = 0.0; % d6
