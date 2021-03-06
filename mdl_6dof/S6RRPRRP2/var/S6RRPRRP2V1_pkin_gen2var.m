% Umwandlung der Kinematikparameter von S6RRPRRP2 zu S6RRPRRP2V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPRRP2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d4 d5 theta3]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPRRP2V1
%   pkin_var=[a3 a4 a6 d1 d4 theta3]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRP2V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 5, 6, 8, 10];
pkin_var = pkin_gen(I_gv);
