% Umwandlung der Kinematikparameter von S5PPRRR1 zu S5PPRRR1V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PPRRR1
%   pkin_gen=[a2 a3 a4 a5 d3 d4 d5 theta1 theta2]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PPRRR1V1
%   pkin_var=[a2 a3 a4 d3 d4 theta1 theta2]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PPRRR1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 5, 6, 8, 9];
pkin_var = pkin_gen(I_gv);
