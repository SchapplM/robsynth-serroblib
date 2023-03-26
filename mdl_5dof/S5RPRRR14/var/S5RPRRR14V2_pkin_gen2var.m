% Umwandlung der Kinematikparameter von S5RPRRR14 zu S5RPRRR14V2
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S5RPRRR14
%   pkin_gen=[a2 a3 a4 a5 alpha2 alpha3 d1 d3 d4 d5 theta2]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RPRRR14V2
%   pkin_var=[alpha2 alpha3 d1 theta2]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RPRRR14V2_pkin_gen2var(pkin_gen)
I_gv = [5, 6, 7, 11];
pkin_var = pkin_gen(I_gv);
