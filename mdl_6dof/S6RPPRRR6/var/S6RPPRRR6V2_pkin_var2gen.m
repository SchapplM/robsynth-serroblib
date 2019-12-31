% Umwandlung der Kinematikparameter von S6RPPRRR6V2 zu S6RPPRRR6
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RPPRRR6V2
%   pkin_var=[a2 a3 a4 a6 d1 d4 d6]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPPRRR6
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d4 d5 d6]
%
% Siehe auch: S6RPPRRR6_structural_kinematic_parameters.m
function pkin_gen = S6RPPRRR6V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 3, 5, 6, 7, 9]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(8) = 0.0; % d5
