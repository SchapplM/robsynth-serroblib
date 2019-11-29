% Umwandlung der Kinematikparameter von S3PPR1V1 zu S3PPR1
% Eingabe:
% pkin_var (2x1) double
%   Kinematikparameter (pkin) von S3PPR1V1
%   pkin_var=[a2 a3]
% Ausgabe:
% pkin_gen (3x1) double
%   Kinematikparameter (pkin) von S3PPR1
%   pkin_gen=[a2 a3 d3]
%
% Siehe auch: S3PPR1_structural_kinematic_parameters.m
function pkin_gen = S3PPR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(3,1);
pkin_gen([1, 2]) = pkin_var;

pkin_gen(3) = 0.0; % d3
