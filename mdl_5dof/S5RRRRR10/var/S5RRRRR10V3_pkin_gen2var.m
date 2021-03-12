% Umwandlung der Kinematikparameter von S5RRRRR10 zu S5RRRRR10V3
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5RRRRR10
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S5RRRRR10V3
%   pkin_var=[a4 d1 d4]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRR10V3_pkin_gen2var(pkin_gen)
I_gv = [3, 6, 9];
pkin_var = pkin_gen(I_gv);
