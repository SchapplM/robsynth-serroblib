% Umwandlung der Kinematikparameter von S4RRPP4 zu S4RRPP4V1
% Eingabe:
% pkin_gen (5x1) double
%   Kinematikparameter (pkin) von S4RRPP4
%   pkin_gen=[a2 a3 a4 d1 d2]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RRPP4V1
%   pkin_var=[a3 a4 d1]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRPP4V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4];
pkin_var = pkin_gen(I_gv);
