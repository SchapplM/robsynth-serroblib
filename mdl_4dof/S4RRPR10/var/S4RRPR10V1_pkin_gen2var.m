% Umwandlung der Kinematikparameter von S4RRPR10 zu S4RRPR10V1
% Eingabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RRPR10
%   pkin_gen=[a2 a3 a4 d1 d2 d4]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4RRPR10V1
%   pkin_var=[a3 a4 d1 d4]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRPR10V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 6];
pkin_var = pkin_gen(I_gv);
