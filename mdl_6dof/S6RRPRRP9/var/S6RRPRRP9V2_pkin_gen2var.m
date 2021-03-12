% Umwandlung der Kinematikparameter von S6RRPRRP9 zu S6RRPRRP9V2
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRP9
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 theta3]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPRRP9V2
%   pkin_var=[a3 a4 a6 d1 d4 theta3]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRP9V2_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 5, 7, 9, 11];
pkin_var = pkin_gen(I_gv);
