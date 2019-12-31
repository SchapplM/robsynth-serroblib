% Umwandlung der Kinematikparameter von S6RRRRPP3V2 zu S6RRRRPP3
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRPP3V2
%   pkin_var=[a2 a3 a5 a6 d1 d2 d3]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RRRRPP3
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4]
%
% Siehe auch: S6RRRRPP3_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPP3V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 4, 5, 6, 7, 8]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(9) = 0.0; % d4
