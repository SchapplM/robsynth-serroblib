% Umwandlung der Kinematikparameter von S6RRRPRR15 zu S6RRRPRR15V8
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRPRR15
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d5 d6]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S6RRRPRR15V8
%   pkin_var=[a3 alpha3 d1 d3]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRR15V8_pkin_gen2var(pkin_gen)
I_gv = [2, 7, 8, 10];
pkin_var = pkin_gen(I_gv);
