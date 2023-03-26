% Umwandlung der Kinematikparameter von S4RRRP7 zu S4RRRP7V3
% Eingabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RRRP7
%   pkin_gen=[a2 a3 a4 d1 d2 d3]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4RRRP7V3
%   pkin_var=[a2 a3 d1 d2 d3]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRRP7V3_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 5, 6];
pkin_var = pkin_gen(I_gv);
