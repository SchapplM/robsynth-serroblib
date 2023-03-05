% Umwandlung der Kinematikparameter von S5RRRRP9 zu S5RRRRP9V2
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRRP9
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S5RRRRP9V2
%   pkin_var=[a4 d1 d4]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRP9V2_pkin_gen2var(pkin_gen)
I_gv = [3, 5, 8];
pkin_var = pkin_gen(I_gv);
