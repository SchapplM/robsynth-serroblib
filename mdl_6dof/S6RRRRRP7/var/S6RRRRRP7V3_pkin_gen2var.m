% Umwandlung der Kinematikparameter von S6RRRRRP7 zu S6RRRRRP7V3
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRRRP7
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S6RRRRRP7V3
%   pkin_var=[a4 a6 d1 d4]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRP7V3_pkin_gen2var(pkin_gen)
I_gv = [3, 5, 7, 10];
pkin_var = pkin_gen(I_gv);
