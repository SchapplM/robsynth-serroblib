% Umwandlung der Kinematikparameter von S5RRPRP9 zu S5RRPRP9V1
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRPRP9
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4 theta3]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRPRP9V1
%   pkin_var=[a3 a4 a5 d1 d4 theta3]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPRP9V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 7, 8];
pkin_var = pkin_gen(I_gv);
