% Umwandlung der Kinematikparameter von S3PRR1V1 zu S3PRR1
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S3PRR1V1
% Ausgabe:
% pkin_gen (4x1) double
%   Kinematikparameter (pkin) von S3PRR1
function pkin_gen = S3PRR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(4,1);
pkin_gen([1, 2, 3]) = pkin_var;
pkin_gen(4) = 0.0;
