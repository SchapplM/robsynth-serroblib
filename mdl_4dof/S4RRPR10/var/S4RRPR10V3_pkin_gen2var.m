% Umwandlung der Kinematikparameter von S4RRPR10 zu S4RRPR10V3
% Eingabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RRPR10
%   pkin_gen=[a2 a3 a4 d1 d2 d4]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RRPR10V3
%   pkin_var=[a2 d1 d2]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRPR10V3_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5];
pkin_var = pkin_gen(I_gv);
