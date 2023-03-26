% Umwandlung der Kinematikparameter von S4RPRR3 zu S4RPRR3V1
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RPRR3
%   pkin_gen=[a2 a3 a4 d1 d3 d4 theta2]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4RPRR3V1
%   pkin_var=[a4 d1 d4 theta2]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RPRR3V1_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 6, 7];
pkin_var = pkin_gen(I_gv);
