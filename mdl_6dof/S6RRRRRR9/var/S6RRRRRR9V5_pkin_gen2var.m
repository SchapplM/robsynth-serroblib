% Umwandlung der Kinematikparameter von S6RRRRRR9 zu S6RRRRRR9V5
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRRR9
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d5 d6]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RRRRRR9V5
%   pkin_var=[a4 a6 d1 d4 d6]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRR9V5_pkin_gen2var(pkin_gen)
I_gv = [3, 5, 8, 11, 13];
pkin_var = pkin_gen(I_gv);
