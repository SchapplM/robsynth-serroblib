% Umwandlung der Kinematikparameter von S5PRRRR10 zu S5PRRRR10V1
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S5PRRRR10
%   pkin_gen=[a2 a3 a4 a5 alpha2 alpha3 d2 d3 d4 d5 theta1]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PRRRR10V1
%   pkin_var=[a2 a3 alpha2 alpha3 d2 d3 theta1]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PRRRR10V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 5, 6, 7, 8, 11];
pkin_var = pkin_gen(I_gv);
