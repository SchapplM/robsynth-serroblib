% Umwandlung der Kinematikparameter von S4RRPR3 zu S4RRPR3V1
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RRPR3
%   pkin_gen=[a2 a3 a4 d1 d2 d4 theta3]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4RRPR3V1
%   pkin_var=[a2 d1 d2 theta3]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRPR3V1_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 7];
pkin_var = pkin_gen(I_gv);
