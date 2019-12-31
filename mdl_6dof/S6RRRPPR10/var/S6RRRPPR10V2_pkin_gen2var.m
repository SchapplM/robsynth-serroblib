% Umwandlung der Kinematikparameter von S6RRRPPR10 zu S6RRRPPR10V2
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRPPR10
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d6 theta5]
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRRPPR10V2
%   pkin_var=[a2 a4 a5 a6 alpha2 d1 d2 d6 theta5]
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPPR10V2_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4, 5, 6, 7, 8, 10, 11];
pkin_var = pkin_gen(I_gv);
