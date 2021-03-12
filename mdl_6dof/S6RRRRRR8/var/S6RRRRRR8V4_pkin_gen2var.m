% Umwandlung der Kinematikparameter von S6RRRRRR8 zu S6RRRRRR8V4
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRRR8
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d5 d6]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRRR8V4
%   pkin_var=[a2 a5 alpha2 d1 d2 d5]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRR8V4_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 6, 8, 9, 12];
pkin_var = pkin_gen(I_gv);
