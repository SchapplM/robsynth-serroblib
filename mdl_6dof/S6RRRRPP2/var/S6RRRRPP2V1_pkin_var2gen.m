% Umwandlung der Kinematikparameter von S6RRRRPP2V1 zu S6RRRRPP2
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RRRRPP2V1
%   pkin_var=[a3 a5 a6 d1 d3]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RRRRPP2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4]
%
% Siehe auch: S6RRRRPP2_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPP2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([2, 4, 5, 6, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(3) = 0.0; % a4
pkin_gen(7) = 0.0; % d2
pkin_gen(9) = 0.0; % d4
