% Umwandlung der Kinematikparameter von S5PRPRR2 zu S5PRPRR2V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PRPRR2
%   pkin_gen=[a2 a3 a4 a5 d2 d4 d5 theta1 theta3]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PRPRR2V1
%   pkin_var=[a2 a3 a4 d2 d4 theta1 theta3]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PRPRR2V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 5, 6, 8, 9];
pkin_var = pkin_gen(I_gv);
