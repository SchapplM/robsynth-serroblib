% Umwandlung der Kinematikparameter von S6PRPRRP3 zu S6PRPRRP3V1
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6PRPRRP3
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d4 d5 theta1 theta3]
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PRPRRP3V1
%   pkin_var=[a2 a3 a4 a6 alpha2 d2 d4 theta1 theta3]
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRPRRP3V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 5, 6, 7, 8, 10, 11];
pkin_var = pkin_gen(I_gv);
