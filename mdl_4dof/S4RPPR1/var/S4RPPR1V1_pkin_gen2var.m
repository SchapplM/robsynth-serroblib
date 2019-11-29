% Umwandlung der Kinematikparameter von S4RPPR1 zu S4RPPR1V1
% Eingabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RPPR1
%   pkin_gen=[a2 a3 a4 d1 d4 theta2]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4RPPR1V1
%   pkin_var=[a2 a3 a4 d1 theta2]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RPPR1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 6];
pkin_var = pkin_gen(I_gv);
