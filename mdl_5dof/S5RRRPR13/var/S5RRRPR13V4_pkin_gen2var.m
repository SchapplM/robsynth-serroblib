% Umwandlung der Kinematikparameter von S5RRRPR13 zu S5RRRPR13V4
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRPR13
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d5]
% Ausgabe:
% pkin_var (1x1) double
%   Kinematikparameter (pkin) von S5RRRPR13V4
%   pkin_var=[d1]
% I_gv (1x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRPR13V4_pkin_gen2var(pkin_gen)
I_gv = [6];
pkin_var = pkin_gen(I_gv);
