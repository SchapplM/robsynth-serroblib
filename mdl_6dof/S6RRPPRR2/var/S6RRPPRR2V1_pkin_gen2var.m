% Umwandlung der Kinematikparameter von S6RRPPRR2 zu S6RRPPRR2V1
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPPRR2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d5 d6 theta3 theta4]
% Ausgabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6RRPPRR2V1
%   pkin_var=[a2 a3 a4 a5 a6 d1 d2 d5 theta3 theta4]
% I_gv (10x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPPRR2V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11];
pkin_var = pkin_gen(I_gv);
