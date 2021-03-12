% Umwandlung der Kinematikparameter von S6RRRPRR15 zu S6RRRPRR15V3
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRPRR15
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d5 d6]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRPRR15V3
%   pkin_var=[a2 a4 a5 alpha2 d1 d2 d5]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRR15V3_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4, 6, 8, 9, 11];
pkin_var = pkin_gen(I_gv);
