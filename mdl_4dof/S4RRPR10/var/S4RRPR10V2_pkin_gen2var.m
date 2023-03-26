% Umwandlung der Kinematikparameter von S4RRPR10 zu S4RRPR10V2
% Eingabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RRPR10
%   pkin_gen=[a2 a3 a4 d1 d2 d4]
% Ausgabe:
% pkin_var (1x1) double
%   Kinematikparameter (pkin) von S4RRPR10V2
%   pkin_var=[d1]
% I_gv (1x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRPR10V2_pkin_gen2var(pkin_gen)
I_gv = [4];
pkin_var = pkin_gen(I_gv);
