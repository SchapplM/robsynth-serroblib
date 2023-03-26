% Umwandlung der Kinematikparameter von S5RPRRR2 zu S5RPRRR2V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RPRRR2
%   pkin_gen=[a2 a3 a4 a5 d1 d3 d4 d5 theta2]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RPRRR2V1
%   pkin_var=[a4 a5 d1 d4 d5 theta2]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RPRRR2V1_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 5, 7, 8, 9];
pkin_var = pkin_gen(I_gv);
