% Umwandlung der Kinematikparameter von S6RRRRPP2 zu S6RRRRPP2V2
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RRRRPP2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRPP2V2
%   pkin_var=[a2 a3 a5 a6 d1 d2 d3]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPP2V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 5, 6, 7, 8];
pkin_var = pkin_gen(I_gv);
