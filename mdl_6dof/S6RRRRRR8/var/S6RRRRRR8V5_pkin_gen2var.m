% Umwandlung der Kinematikparameter von S6RRRRRR8 zu S6RRRRRR8V5
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRRR8
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d5 d6]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRRR8V5
%   pkin_var=[a3 a5 alpha3 d1 d3 d5]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRR8V5_pkin_gen2var(pkin_gen)
I_gv = [2, 4, 7, 8, 10, 12];
pkin_var = pkin_gen(I_gv);
