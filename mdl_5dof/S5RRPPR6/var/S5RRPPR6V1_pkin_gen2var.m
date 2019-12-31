% Umwandlung der Kinematikparameter von S5RRPPR6 zu S5RRPPR6V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPPR6
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d5 theta3 theta4]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRPPR6V1
%   pkin_var=[a3 a4 a5 d1 d5 theta3 theta4]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPPR6V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 7, 8, 9];
pkin_var = pkin_gen(I_gv);
