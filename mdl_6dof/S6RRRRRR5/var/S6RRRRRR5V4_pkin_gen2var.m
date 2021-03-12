% Umwandlung der Kinematikparameter von S6RRRRRR5 zu S6RRRRRR5V4
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRRR5
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d4 d5 d6]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RRRRRR5V4
%   pkin_var=[a4 a5 d1 d4 d5]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRR5V4_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 7, 10, 11];
pkin_var = pkin_gen(I_gv);
