% Umwandlung der Kinematikparameter von S5PRRPR7 zu S5PRRPR7V4
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5PRRPR7
%   pkin_gen=[a2 a3 a4 a5 alpha2 d2 d3 d5 theta1 theta4]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5PRRPR7V4
%   pkin_var=[a2 a4 a5 d2 d5 theta1]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PRRPR7V4_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4, 6, 8, 9];
pkin_var = pkin_gen(I_gv);
