% Umwandlung der Kinematikparameter von S5PRPRR8 zu S5PRPRR8V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PRPRR8
%   pkin_gen=[a2 a3 a4 a5 alpha2 d2 d4 d5 theta1]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PRPRR8V1
%   pkin_var=[a2 a3 a4 alpha2 d2 d4 theta1]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PRPRR8V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 5, 6, 7, 9];
pkin_var = pkin_gen(I_gv);
