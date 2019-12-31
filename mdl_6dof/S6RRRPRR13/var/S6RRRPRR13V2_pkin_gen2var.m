% Umwandlung der Kinematikparameter von S6RRRPRR13 zu S6RRRPRR13V2
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRPRR13
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d5 d6 theta4]
% Ausgabe:
% pkin_var (11x1) double
%   Kinematikparameter (pkin) von S6RRRPRR13V2
%   pkin_var=[a2 a3 a4 a5 alpha2 alpha3 d1 d2 d3 d5 theta4]
% I_gv (11x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRR13V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 13];
pkin_var = pkin_gen(I_gv);
