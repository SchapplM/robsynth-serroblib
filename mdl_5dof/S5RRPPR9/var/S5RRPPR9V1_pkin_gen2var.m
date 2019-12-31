% Umwandlung der Kinematikparameter von S5RRPPR9 zu S5RRPPR9V1
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S5RRPPR9
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d5]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRPPR9V1
%   pkin_var=[a3 a4 a5 d1 d5]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPPR9V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 7];
pkin_var = pkin_gen(I_gv);
