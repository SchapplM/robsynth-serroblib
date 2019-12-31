% Umwandlung der Kinematikparameter von S4RRPR6 zu S4RRPR6V1
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RRPR6
%   pkin_gen=[a2 a3 a4 d1 d2 d4 theta3]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4RRPR6V1
%   pkin_var=[a3 a4 d1 d4 theta3]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRPR6V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 6, 7];
pkin_var = pkin_gen(I_gv);
