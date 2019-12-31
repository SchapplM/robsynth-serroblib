% Umwandlung der Kinematikparameter von S5RPRRR14 zu S5RPRRR14V1
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S5RPRRR14
%   pkin_gen=[a2 a3 a4 a5 alpha2 alpha3 d1 d3 d4 d5 theta2]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RPRRR14V1
%   pkin_var=[a2 a3 alpha2 alpha3 d1 d3 theta2]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RPRRR14V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 5, 6, 7, 8, 11];
pkin_var = pkin_gen(I_gv);
