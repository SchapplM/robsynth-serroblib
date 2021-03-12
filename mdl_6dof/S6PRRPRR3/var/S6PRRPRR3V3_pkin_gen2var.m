% Umwandlung der Kinematikparameter von S6PRRPRR3 zu S6PRRPRR3V3
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6PRRPRR3
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d2 d3 d5 d6 theta1 theta4]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6PRRPRR3V3
%   pkin_var=[a2 a4 a5 alpha2 d2 d5 theta1 theta4]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRRPRR3V3_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4, 6, 8, 10, 12, 13];
pkin_var = pkin_gen(I_gv);
