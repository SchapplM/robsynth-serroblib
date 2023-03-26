% Umwandlung der Kinematikparameter von S6RRPRRR14 zu S6RRPRRR14V5
% Eingabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RRPRRR14
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d1 d2 d4 d5 d6 theta3]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRPRRR14V5
%   pkin_var=[a2 alpha2 alpha3 alpha4 d1 d2 theta3]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR14V5_pkin_gen2var(pkin_gen)
I_gv = [1, 6, 7, 8, 9, 10, 14];
pkin_var = pkin_gen(I_gv);
