% Umwandlung der Kinematikparameter von S6RRRPRR14 zu S6RRRPRR14V3
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRPRR14
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d5 d6]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRPRR14V3
%   pkin_var=[a2 a6 alpha2 d1 d2 d6]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRR14V3_pkin_gen2var(pkin_gen)
I_gv = [1, 5, 6, 7, 8, 11];
pkin_var = pkin_gen(I_gv);
