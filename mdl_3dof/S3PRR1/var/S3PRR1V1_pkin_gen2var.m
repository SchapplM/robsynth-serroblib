% Umwandlung der Kinematikparameter von S3PRR1 zu S3PRR1V1
% Eingabe:
% pkin_gen (4x1) double
%   Kinematikparameter (pkin) von S3PRR1
%   pkin_gen=[a2 a3 d2 d3]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S3PRR1V1
%   pkin_var=[a2 a3 d2]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S3PRR1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3];
pkin_var = pkin_gen(I_gv);
