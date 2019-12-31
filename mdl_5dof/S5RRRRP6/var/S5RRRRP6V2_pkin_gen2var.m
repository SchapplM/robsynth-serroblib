% Umwandlung der Kinematikparameter von S5RRRRP6 zu S5RRRRP6V2
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRRP6
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRRRP6V2
%   pkin_var=[a2 a3 a5 d1 d2 d3]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRP6V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 5, 6, 7];
pkin_var = pkin_gen(I_gv);
