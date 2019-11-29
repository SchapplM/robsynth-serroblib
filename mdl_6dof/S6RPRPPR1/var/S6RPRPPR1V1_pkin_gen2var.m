% Umwandlung der Kinematikparameter von S6RPRPPR1 zu S6RPRPPR1V1
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RPRPPR1
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d6 theta2 theta4 theta5]
% Ausgabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6RPRPPR1V1
%   pkin_var=[a2 a3 a4 a5 a6 d1 d3 theta2 theta4 theta5]
% I_gv (10x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRPPR1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 9, 10, 11];
pkin_var = pkin_gen(I_gv);
