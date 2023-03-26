% Umwandlung der Kinematikparameter von S4RRPR9 zu S4RRPR9V2
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RRPR9
%   pkin_gen=[a2 a3 a4 d1 d2 d4 theta3]
% Ausgabe:
% pkin_var (2x1) double
%   Kinematikparameter (pkin) von S4RRPR9V2
%   pkin_var=[d1 theta3]
% I_gv (2x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRPR9V2_pkin_gen2var(pkin_gen)
I_gv = [4, 7];
pkin_var = pkin_gen(I_gv);
