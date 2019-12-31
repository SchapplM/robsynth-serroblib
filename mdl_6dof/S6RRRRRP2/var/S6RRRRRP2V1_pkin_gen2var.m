% Umwandlung der Kinematikparameter von S6RRRRRP2 zu S6RRRRRP2V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRRRP2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRRP2V1
%   pkin_var=[a3 a4 a6 d1 d3 d4]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRP2V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 5, 6, 8, 9];
pkin_var = pkin_gen(I_gv);
