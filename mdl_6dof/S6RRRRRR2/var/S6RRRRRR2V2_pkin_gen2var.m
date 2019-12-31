% Umwandlung der Kinematikparameter von S6RRRRRR2 zu S6RRRRRR2V2
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRRRR2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d5 d6]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRRR2V2
%   pkin_var=[a3 a4 a6 d1 d3 d4 d6]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRR2V2_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 5, 6, 8, 9, 11];
pkin_var = pkin_gen(I_gv);
