% Umwandlung der Kinematikparameter von S3PPR1 zu S3PPR1V1
% Eingabe:
% pkin_gen (3x1) double
%   Kinematikparameter (pkin) von S3PPR1
% Ausgabe:
% pkin_var (2x1) double
%   Kinematikparameter (pkin) von S3PPR1V1
% I_gv (2x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S3PPR1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2];
pkin_var = pkin_gen(I_gv);
