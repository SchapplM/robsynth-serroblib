% Umwandlung der Kinematikparameter von S6RRRPRR9 zu S6RRRPRR9V4
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRPRR9
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d5 d6 theta4]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRPRR9V4
%   pkin_var=[a3 a4 a5 alpha3 d1 d3 d5 theta4]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRR9V4_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 7, 8, 10, 11, 13];
pkin_var = pkin_gen(I_gv);
