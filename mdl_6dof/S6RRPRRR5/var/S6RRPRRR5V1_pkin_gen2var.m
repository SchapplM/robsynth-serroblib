% Umwandlung der Kinematikparameter von S6RRPRRR5 zu S6RRPRRR5V1
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRPRRR5
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 d6 theta3]
% Ausgabe:
% pkin_var (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR5V1
%   pkin_var=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 theta3]
% I_gv (11x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR5V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12];
pkin_var = pkin_gen(I_gv);
