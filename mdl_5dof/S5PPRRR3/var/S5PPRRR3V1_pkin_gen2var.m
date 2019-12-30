% Umwandlung der Kinematikparameter von S5PPRRR3 zu S5PPRRR3V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PPRRR3
%   pkin_gen=[a2 a3 a4 a5 d3 d4 d5 theta1 theta2]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S5PPRRR3V1
%   pkin_var=[a2 a3 a4 a5 d3 d4 d5 theta1]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PPRRR3V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 8];
pkin_var = pkin_gen(I_gv);
