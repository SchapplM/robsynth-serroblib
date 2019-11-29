% Umwandlung der Kinematikparameter von S4RPRR1 zu S4RPRR1V1
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RPRR1
%   pkin_gen=[a2 a3 a4 d1 d3 d4 theta2]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S4RPRR1V1
%   pkin_var=[a2 a3 a4 d1 d3 theta2]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RPRR1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 7];
pkin_var = pkin_gen(I_gv);
