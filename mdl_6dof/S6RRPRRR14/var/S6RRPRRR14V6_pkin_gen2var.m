% Umwandlung der Kinematikparameter von S6RRPRRR14 zu S6RRPRRR14V6
% Eingabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RRPRRR14
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d1 d2 d4 d5 d6 theta3]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S6RRPRRR14V6
%   pkin_var=[alpha3 alpha4 d1 theta3]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR14V6_pkin_gen2var(pkin_gen)
I_gv = [7, 8, 9, 14];
pkin_var = pkin_gen(I_gv);
