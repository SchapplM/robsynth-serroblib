% Umwandlung der Kinematikparameter von S6RPPRRR5V1 zu S6RPPRRR5
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RPPRRR5V1
%   pkin_var=[a2 a3 a4 a5 a6 d1 d4 d5]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPPRRR5
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d4 d5 d6]
%
% Siehe auch: S6RPPRRR5_structural_kinematic_parameters.m
function pkin_gen = S6RPPRRR5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8]) = pkin_var;

pkin_gen(9) = 0.0; % d6
