% Umwandlung der Kinematikparameter von S6RRPRRR13 zu S6RRPRRR13V6
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR13
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 d6]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRPRRR13V6
%   pkin_var=[a2 a5 a6 alpha2 d1 d2 d5 d6]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR13V6_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 7, 8, 10, 11];
pkin_var = pkin_gen(I_gv);
