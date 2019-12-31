% Umwandlung der Kinematikparameter von S6RPPPRR4 zu S6RPPPRR4V2
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPPPRR4
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d5 d6 theta3]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RPPPRR4V2
%   pkin_var=[a2 a3 a4 a5 d1 d5 theta3]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPPPRR4V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 6, 7, 9];
pkin_var = pkin_gen(I_gv);
