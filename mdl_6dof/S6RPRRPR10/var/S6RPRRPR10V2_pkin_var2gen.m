% Umwandlung der Kinematikparameter von S6RPRRPR10V2 zu S6RPRRPR10
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RPRRPR10V2
%   pkin_var=[a2 a3 a5 a6 d1 d3 d6]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPRRPR10
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 d6]
%
% Siehe auch: S6RPRRPR10_structural_kinematic_parameters.m
function pkin_gen = S6RPRRPR10V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 4, 5, 6, 7, 9]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(8) = 0.0; % d4
