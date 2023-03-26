% Umwandlung der Kinematikparameter von S5RPRRR14 zu S5RPRRR14V3
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S5RPRRR14
%   pkin_gen=[a2 a3 a4 a5 alpha2 alpha3 d1 d3 d4 d5 theta2]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S5RPRRR14V3
%   pkin_var=[a4 a5 alpha2 alpha3 d1 d4 d5 theta2]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RPRRR14V3_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 5, 6, 7, 9, 10, 11];
pkin_var = pkin_gen(I_gv);
