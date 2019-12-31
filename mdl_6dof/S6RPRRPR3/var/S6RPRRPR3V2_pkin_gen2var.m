% Umwandlung der Kinematikparameter von S6RPRRPR3 zu S6RPRRPR3V2
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RPRRPR3
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 d6 theta2]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RPRRPR3V2
%   pkin_var=[a2 a3 a5 a6 d1 d3 d6 theta2]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRPR3V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 5, 6, 7, 9, 10];
pkin_var = pkin_gen(I_gv);
