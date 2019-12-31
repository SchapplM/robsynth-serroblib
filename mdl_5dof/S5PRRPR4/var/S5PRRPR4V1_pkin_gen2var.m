% Umwandlung der Kinematikparameter von S5PRRPR4 zu S5PRRPR4V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PRRPR4
%   pkin_gen=[a2 a3 a4 a5 d2 d3 d5 theta1 theta4]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PRRPR4V1
%   pkin_var=[a2 a4 a5 d2 d5 theta1 theta4]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PRRPR4V1_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4, 5, 7, 8, 9];
pkin_var = pkin_gen(I_gv);
