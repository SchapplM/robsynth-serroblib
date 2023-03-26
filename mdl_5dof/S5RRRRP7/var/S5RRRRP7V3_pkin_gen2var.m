% Umwandlung der Kinematikparameter von S5RRRRP7 zu S5RRRRP7V3
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRRP7
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S5RRRRP7V3
%   pkin_var=[a3 d1 d3]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRP7V3_pkin_gen2var(pkin_gen)
I_gv = [2, 5, 7];
pkin_var = pkin_gen(I_gv);
