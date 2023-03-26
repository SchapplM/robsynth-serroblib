% Umwandlung der Kinematikparameter von S6RRPRRR8 zu S6RRPRRR8V3
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR8
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d4 d5 d6 theta3]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRPRRR8V3
%   pkin_var=[a2 a5 a6 d1 d2 d5 d6 theta3]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR8V3_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 7, 9, 10, 11];
pkin_var = pkin_gen(I_gv);
