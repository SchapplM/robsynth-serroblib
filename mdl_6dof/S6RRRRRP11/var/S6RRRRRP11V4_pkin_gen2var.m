% Umwandlung der Kinematikparameter von S6RRRRRP11 zu S6RRRRRP11V4
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRRP11
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S6RRRRRP11V4
%   pkin_var=[a4 a6 d1 d4]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRP11V4_pkin_gen2var(pkin_gen)
I_gv = [3, 5, 8, 11];
pkin_var = pkin_gen(I_gv);
