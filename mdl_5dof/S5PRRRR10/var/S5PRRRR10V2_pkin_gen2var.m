% Umwandlung der Kinematikparameter von S5PRRRR10 zu S5PRRRR10V2
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S5PRRRR10
%   pkin_gen=[a2 a3 a4 a5 alpha2 alpha3 d2 d3 d4 d5 theta1]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5PRRRR10V2
%   pkin_var=[a2 a5 alpha2 d2 d5 theta1]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PRRRR10V2_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 7, 10, 11];
pkin_var = pkin_gen(I_gv);
