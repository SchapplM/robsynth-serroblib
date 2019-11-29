% Umwandlung der Kinematikparameter von S3RRR1V1 zu S3RRR1
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S3RRR1V1
%   pkin_var=[a2 a3 d1 d2]
% Ausgabe:
% pkin_gen (5x1) double
%   Kinematikparameter (pkin) von S3RRR1
%   pkin_gen=[a2 a3 d1 d2 d3]
%
% Siehe auch: S3RRR1_structural_kinematic_parameters.m
function pkin_gen = S3RRR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(5,1);
pkin_gen([1, 2, 3, 4]) = pkin_var;

pkin_gen(5) = 0.0; % d3
