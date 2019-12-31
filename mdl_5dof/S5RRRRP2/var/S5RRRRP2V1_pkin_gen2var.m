% Umwandlung der Kinematikparameter von S5RRRRP2 zu S5RRRRP2V1
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRRP2
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRRRP2V1
%   pkin_var=[a2 a4 a5 d1 d2 d4]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRP2V1_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4, 5, 6, 8];
pkin_var = pkin_gen(I_gv);
