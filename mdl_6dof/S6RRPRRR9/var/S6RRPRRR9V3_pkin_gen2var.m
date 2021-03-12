% Umwandlung der Kinematikparameter von S6RRPRRR9 zu S6RRPRRR9V3
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRPRRR9
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 d6 theta3]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRPRRR9V3
%   pkin_var=[a3 a4 a5 d1 d4 d5 theta3]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR9V3_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 7, 9, 10, 12];
pkin_var = pkin_gen(I_gv);
