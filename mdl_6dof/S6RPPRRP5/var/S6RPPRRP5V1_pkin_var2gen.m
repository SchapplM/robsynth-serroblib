% Umwandlung der Kinematikparameter von S6RPPRRP5V1 zu S6RPPRRP5
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RPPRRP5V1
%   pkin_var=[a2 a3 a4 a6 d1 d4]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S6RPPRRP5
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d4 d5]
%
% Siehe auch: S6RPPRRP5_structural_kinematic_parameters.m
function pkin_gen = S6RPPRRP5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 2, 3, 5, 6, 7]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(8) = 0.0; % d5
