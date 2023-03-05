% Umwandlung der Kinematikparameter von S5RRPRR14 zu S5RRPRR14V3
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5RRPRR14
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d4 d5 theta3]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRPRR14V3
%   pkin_var=[a2 alpha2 d1 d2 theta3]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPRR14V3_pkin_gen2var(pkin_gen)
I_gv = [1, 5, 6, 7, 10];
pkin_var = pkin_gen(I_gv);
