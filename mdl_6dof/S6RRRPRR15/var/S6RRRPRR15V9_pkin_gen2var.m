% Umwandlung der Kinematikparameter von S6RRRPRR15 zu S6RRRPRR15V9
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRPRR15
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d5 d6]
% Ausgabe:
% pkin_var (1x1) double
%   Kinematikparameter (pkin) von S6RRRPRR15V9
%   pkin_var=[d1]
% I_gv (1x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRR15V9_pkin_gen2var(pkin_gen)
I_gv = [8];
pkin_var = pkin_gen(I_gv);
