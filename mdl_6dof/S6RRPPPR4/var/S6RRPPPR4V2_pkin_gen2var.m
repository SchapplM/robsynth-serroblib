% Umwandlung der Kinematikparameter von S6RRPPPR4 zu S6RRPPPR4V2
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RRPPPR4
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d6 theta4]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRPPPR4V2
%   pkin_var=[a3 a4 a5 a6 d1 d6 theta4]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPPPR4V2_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 6, 8, 9];
pkin_var = pkin_gen(I_gv);
