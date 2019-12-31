% Umwandlung der Kinematikparameter von S5PPPRR2 zu S5PPPRR2V4
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PPPRR2
%   pkin_gen=[a2 a3 a4 a5 d4 d5 theta1 theta2 theta3]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5PPPRR2V4
%   pkin_var=[a2 a3 a4 d4 theta1]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PPPRR2V4_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 5, 7];
pkin_var = pkin_gen(I_gv);
