% Umwandlung der Kinematikparameter von S6RRPRRR7 zu S6RRPRRR7V2
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPRRR7
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d4 d5 d6]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPRRR7V2
%   pkin_var=[a3 a4 a6 d1 d4 d6]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR7V2_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 5, 6, 8, 10];
pkin_var = pkin_gen(I_gv);
