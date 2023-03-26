% Umwandlung der Kinematikparameter von S6RRRPRR4 zu S6RRRPRR4V4
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRPRR4
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d5 d6 theta4]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRPRR4V4
%   pkin_var=[a2 a3 a6 d1 d2 d3 d6 theta4]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRR4V4_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 5, 6, 7, 8, 10, 11];
pkin_var = pkin_gen(I_gv);
