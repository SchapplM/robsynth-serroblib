% Umwandlung der Kinematikparameter von S4PPRR3 zu S4PPRR3V1
% Eingabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4PPRR3
%   pkin_gen=[a2 a3 a4 d3 d4 theta1]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4PPRR3V1
%   pkin_var=[a2 a3 d3 theta1]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4PPRR3V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 6];
pkin_var = pkin_gen(I_gv);
