% Umwandlung der Kinematikparameter von S6RRRRRP2 zu S6RRRRRP2V4
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRRRP2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRRP2V4
%   pkin_var=[a2 a3 a4 d1 d2 d3 d4]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRP2V4_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 8, 9];
pkin_var = pkin_gen(I_gv);
