% Umwandlung der Kinematikparameter von S5RRRRR10 zu S5RRRRR10V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5RRRRR10
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRRRR10V1
%   pkin_var=[a2 a4 alpha2 d1 d2 d4]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRR10V1_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 5, 6, 7, 9];
pkin_var = pkin_gen(I_gv);
