% Umwandlung der Kinematikparameter von S5RRPRR13 zu S5RRPRR13V2
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPRR13
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4 d5 theta3]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRPRR13V2
%   pkin_var=[a5 d1 d5 theta3]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPRR13V2_pkin_gen2var(pkin_gen)
I_gv = [4, 5, 8, 9];
pkin_var = pkin_gen(I_gv);
