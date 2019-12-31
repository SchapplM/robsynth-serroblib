% Umwandlung der Kinematikparameter von S6RRRRPR6V2 zu S6RRRRPR6
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRPR6V2
%   pkin_var=[a4 a5 a6 d1 d4 d6 theta5]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRRPR6
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d6 theta5]
%
% Siehe auch: S6RRRRPR6_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPR6V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([3, 4, 5, 6, 9, 10, 11]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(7) = 0.0; % d2
pkin_gen(8) = 0.0; % d3
