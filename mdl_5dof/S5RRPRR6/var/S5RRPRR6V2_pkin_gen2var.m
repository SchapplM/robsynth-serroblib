% Umwandlung der Kinematikparameter von S5RRPRR6 zu S5RRPRR6V2
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPRR6
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4 d5 theta3]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRPRR6V2
%   pkin_var=[a2 a5 d1 d2 d5]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPRR6V2_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 8];
pkin_var = pkin_gen(I_gv);
