% Umwandlung der Kinematikparameter von S4RRPR7 zu S4RRPR7V3
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RRPR7
%   pkin_gen=[a2 a3 a4 d1 d2 d4 theta3]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4RRPR7V3
%   pkin_var=[a2 d1 d2 theta3]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRPR7V3_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 7];
pkin_var = pkin_gen(I_gv);
