% Umwandlung der Kinematikparameter von S5RRRPR9 zu S5RRRPR9V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRPR9
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d5 theta4]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRRPR9V1
%   pkin_var=[a4 a5 d1 d5 theta4]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRPR9V1_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 5, 8, 9];
pkin_var = pkin_gen(I_gv);
