% Umwandlung der Kinematikparameter von S6RRPRRR13 zu S6RRPRRR13V5
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR13
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 d6]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S6RRPRRR13V5
%   pkin_var=[a6 d1 d6]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR13V5_pkin_gen2var(pkin_gen)
I_gv = [5, 7, 11];
pkin_var = pkin_gen(I_gv);
