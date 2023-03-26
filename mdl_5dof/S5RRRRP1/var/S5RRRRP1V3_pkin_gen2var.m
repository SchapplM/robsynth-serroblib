% Umwandlung der Kinematikparameter von S5RRRRP1 zu S5RRRRP1V3
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRRP1
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRRRP1V3
%   pkin_var=[a2 a3 a4 d1 d2 d3 d4]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRP1V3_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 5, 6, 7, 8];
pkin_var = pkin_gen(I_gv);
