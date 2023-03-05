% Umwandlung der Kinematikparameter von S5RPRRR13 zu S5RPRRR13V3
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RPRRR13
%   pkin_gen=[a2 a3 a4 a5 d1 d3 d4 d5]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RPRRR13V3
%   pkin_var=[a4 a5 d1 d4 d5]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RPRRR13V3_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 5, 7, 8];
pkin_var = pkin_gen(I_gv);
