% Umwandlung der Kinematikparameter von S5RRPRR10 zu S5RRPRR10V5
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5RRPRR10
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d4 d5 theta3]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRPRR10V5
%   pkin_var=[a2 a5 alpha2 d1 d2 d5 theta3]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPRR10V5_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 7, 9, 10];
pkin_var = pkin_gen(I_gv);
