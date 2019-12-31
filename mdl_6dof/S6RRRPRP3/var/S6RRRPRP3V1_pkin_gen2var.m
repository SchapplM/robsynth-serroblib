% Umwandlung der Kinematikparameter von S6RRRPRP3 zu S6RRRPRP3V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRPRP3
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d5 theta4]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRPRP3V1
%   pkin_var=[a3 a4 a5 a6 d1 d3 d5 theta4]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRP3V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 6, 8, 9, 10];
pkin_var = pkin_gen(I_gv);
