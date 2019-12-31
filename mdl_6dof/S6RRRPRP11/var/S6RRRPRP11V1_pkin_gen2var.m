% Umwandlung der Kinematikparameter von S6RRRPRP11 zu S6RRRPRP11V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRPRP11
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d5]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRPRP11V1
%   pkin_var=[a2 a4 a5 a6 alpha2 d1 d2 d5]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRP11V1_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4, 5, 6, 7, 8, 10];
pkin_var = pkin_gen(I_gv);
