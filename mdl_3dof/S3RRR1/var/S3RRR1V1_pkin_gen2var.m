% Umwandlung der Kinematikparameter von S3RRR1 zu S3RRR1V1
% Eingabe:
% pkin_gen (5x1) double
%   Kinematikparameter (pkin) von S3RRR1
%   pkin_gen=[a2 a3 d1 d2 d3]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S3RRR1V1
%   pkin_var=[a2 a3 d1 d2]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S3RRR1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4];
pkin_var = pkin_gen(I_gv);
