% Umwandlung der Kinematikparameter von S6RPRRPP5 zu S6RPRRPP5V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPRRPP5
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 theta2]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RPRRPP5V1
%   pkin_var=[a2 a3 a5 a6 d1 d3 theta2]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRPP5V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 5, 6, 7, 9];
pkin_var = pkin_gen(I_gv);
