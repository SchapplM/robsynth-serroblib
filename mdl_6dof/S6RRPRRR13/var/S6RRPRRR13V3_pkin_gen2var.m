% Umwandlung der Kinematikparameter von S6RRPRRR13 zu S6RRPRRR13V3
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR13
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 d6]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPRRR13V3
%   pkin_var=[a3 a4 a6 d1 d4 d6]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR13V3_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 5, 7, 9, 11];
pkin_var = pkin_gen(I_gv);
