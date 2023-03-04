% Umwandlung der Kinematikparameter von S6RRRRPR7V4 zu S6RRRRPR7
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRPR7V4
%   pkin_var=[a2 a4 alpha2 d1 d2 d4 theta5]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRPR7
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d4 d6 theta5]
%
% Siehe auch: S6RRRRPR7_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPR7V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 3, 6, 7, 8, 10, 12]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(9) = 0.0; % d3
pkin_gen(11) = 0.0; % d6
