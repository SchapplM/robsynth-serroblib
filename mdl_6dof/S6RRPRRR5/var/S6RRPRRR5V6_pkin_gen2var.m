% Umwandlung der Kinematikparameter von S6RRPRRR5 zu S6RRPRRR5V6
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRPRRR5
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 d6 theta3]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S6RRPRRR5V6
%   pkin_var=[a6 d1 d6 theta3]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR5V6_pkin_gen2var(pkin_gen)
I_gv = [5, 7, 11, 12];
pkin_var = pkin_gen(I_gv);
