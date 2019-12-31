% Umwandlung der Kinematikparameter von S6RPRRPP7 zu S6RPRRPP7V1
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S6RPRRPP7
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RPRRPP7V1
%   pkin_var=[a2 a3 a5 a6 d1 d3]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRPP7V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 5, 6, 7];
pkin_var = pkin_gen(I_gv);
