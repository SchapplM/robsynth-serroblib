% Umwandlung der Kinematikparameter von S6RRRRRR10 zu S6RRRRRR10V5
% Eingabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RRRRRR10
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d1 d2 d3 d4 d5 d6]
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRRRRR10V5
%   pkin_var=[a2 a3 a6 alpha2 alpha3 d1 d2 d3 d6]
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRR10V5_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 5, 6, 7, 9, 10, 11, 14];
pkin_var = pkin_gen(I_gv);
