% Umwandlung der Kinematikparameter von S6RRRPRR11V8 zu S6RRRPRR11
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRPRR11V8
%   pkin_var=[a2 a3 a6 alpha2 d1 d2 d3 d6]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRPRR11
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d5 d6]
%
% Siehe auch: S6RRRPRR11_structural_kinematic_parameters.m
function pkin_gen = S6RRRPRR11V8_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 5, 6, 7, 8, 9, 11]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(10) = 0.0; % d5
