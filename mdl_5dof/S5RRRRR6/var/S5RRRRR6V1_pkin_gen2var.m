% Umwandlung der Kinematikparameter von S5RRRRR6 zu S5RRRRR6V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRRR6
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRRRR6V1
%   pkin_var=[a2 a4 a5 d1 d2 d4 d5]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRR6V1_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4, 5, 6, 8, 9];
pkin_var = pkin_gen(I_gv);
