% Umwandlung der Kinematikparameter von S6PRRPPR1 zu S6PRRPPR1V1
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRRPPR1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d3 d6 theta1 theta4 theta5]
% Ausgabe:
% pkin_var (11x1) double
%   Kinematikparameter (pkin) von S6PRRPPR1V1
%   pkin_var=[a2 a3 a4 a5 a6 alpha2 d2 d3 theta1 theta4 theta5]
% I_gv (11x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRRPPR1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12];
pkin_var = pkin_gen(I_gv);
