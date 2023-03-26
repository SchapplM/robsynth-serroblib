% Umwandlung der Kinematikparameter von S6RRRRRP8 zu S6RRRRRP8V7
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRRRP8
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6RRRRRP8V7
%   pkin_var=[a2 a3 a4 a5 alpha2 d1 d2 d3 d4 d5]
% I_gv (10x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRP8V7_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 6, 7, 8, 9, 10, 11];
pkin_var = pkin_gen(I_gv);
