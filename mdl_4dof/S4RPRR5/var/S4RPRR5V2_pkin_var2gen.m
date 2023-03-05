% Umwandlung der Kinematikparameter von S4RPRR5V2 zu S4RPRR5
% Eingabe:
% pkin_var (1x1) double
%   Kinematikparameter (pkin) von S4RPRR5V2
%   pkin_var=[d1]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RPRR5
%   pkin_gen=[a2 a3 a4 d1 d3 d4]
%
% Siehe auch: S4RPRR5_structural_kinematic_parameters.m
function pkin_gen = S4RPRR5V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([4]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(5) = 0.0; % d3
pkin_gen(6) = 0.0; % d4
