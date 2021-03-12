% Umwandlung der Kinematikparameter von S6RRRPRR11 zu S6RRRPRR11V4
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRPRR11
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d5 d6]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S6RRRPRR11V4
%   pkin_var=[a4 a5 d1 d5]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRR11V4_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 7, 10];
pkin_var = pkin_gen(I_gv);
