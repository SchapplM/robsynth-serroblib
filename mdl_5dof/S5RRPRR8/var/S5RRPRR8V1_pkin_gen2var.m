% Umwandlung der Kinematikparameter von S5RRPRR8 zu S5RRPRR8V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPRR8
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4 d5 theta3]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRPRR8V1
%   pkin_var=[a3 a4 d1 d4 theta3]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPRR8V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 5, 7, 9];
pkin_var = pkin_gen(I_gv);
