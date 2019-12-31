% Umwandlung der Kinematikparameter von S5RPPRR11 zu S5RPPRR11V1
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S5RPPRR11
%   pkin_gen=[a2 a3 a4 a5 d1 d4 d5]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RPPRR11V1
%   pkin_var=[a2 a3 a4 d1 d4]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RPPRR11V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 5, 6];
pkin_var = pkin_gen(I_gv);
