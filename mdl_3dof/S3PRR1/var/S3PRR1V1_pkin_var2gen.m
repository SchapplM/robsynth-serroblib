% Umwandlung der Kinematikparameter von S3PRR1V1 zu S3PRR1
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S3PRR1V1
%   pkin_var=[a2 a3 d2]
% Ausgabe:
% pkin_gen (4x1) double
%   Kinematikparameter (pkin) von S3PRR1
%   pkin_gen=[a2 a3 d2 d3]
%
% Siehe auch: S3PRR1_structural_kinematic_parameters.m
function pkin_gen = S3PRR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(4,1);
pkin_gen([1, 2, 3]) = pkin_var;

pkin_gen(4) = 0.0; % d3
