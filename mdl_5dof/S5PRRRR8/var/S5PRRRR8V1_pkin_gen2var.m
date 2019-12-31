% Umwandlung der Kinematikparameter von S5PRRRR8 zu S5PRRRR8V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5PRRRR8
%   pkin_gen=[a2 a3 a4 a5 alpha2 d2 d3 d4 d5 theta1]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5PRRRR8V1
%   pkin_var=[a2 a4 alpha2 d2 d4 theta1]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PRRRR8V1_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 5, 6, 8, 10];
pkin_var = pkin_gen(I_gv);
