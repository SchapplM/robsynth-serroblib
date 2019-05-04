% Umwandlung der Kinematikparameter von S4PPPR3 zu S4PPPR3V1
% Eingabe:
% pkin_gen (5x1) double
%   Kinematikparameter (pkin) von S4PPPR3
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4PPPR3V1
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4PPPR3V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 5];
pkin_var = pkin_gen(I_gv);
