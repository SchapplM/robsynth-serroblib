% Umwandlung der Kinematikparameter von S6PRPRPR7 zu S6PRPRPR7V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6PRPRPR7
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d4 d6 theta1]
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PRPRPR7V1
%   pkin_var=[a2 a3 a4 a5 a6 alpha2 d2 d4 theta1]
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRPRPR7V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 8, 10];
pkin_var = pkin_gen(I_gv);
