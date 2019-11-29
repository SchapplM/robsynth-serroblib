% Umwandlung der Kinematikparameter von S6RPRRRR2 zu S6RPRRRR2V1
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RPRRRR2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 d5 d6 theta2]
% Ausgabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6RPRRRR2V1
%   pkin_var=[a2 a3 a4 a5 a6 d1 d3 d4 d5 theta2]
% I_gv (10x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRRR2V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11];
pkin_var = pkin_gen(I_gv);
