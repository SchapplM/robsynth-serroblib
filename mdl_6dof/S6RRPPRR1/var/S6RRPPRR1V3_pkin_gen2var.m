% Umwandlung der Kinematikparameter von S6RRPPRR1 zu S6RRPPRR1V3
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPPRR1
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d5 d6 theta3]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRPPRR1V3
%   pkin_var=[a2 a3 a4 a5 d1 d2 d5 theta3]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPPRR1V3_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 6, 7, 8, 10];
pkin_var = pkin_gen(I_gv);
