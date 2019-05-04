% Umwandlung der Kinematikparameter von S6RRRRRR10V2 zu S6RRRRRR10
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRRR10V2
% Ausgabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RRRRRR10
function pkin_gen = S6RRRRRR10V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(14,1);
pkin_gen([1, 2, 3, 9, 12, 14]) = pkin_var;
pkin_gen(4) = 0.0;
pkin_gen(5) = 0.0;
pkin_gen(6) = pi/2;
pkin_gen(7) = 0.0;
pkin_gen(8) = pi/2;
pkin_gen(10) = 0.0;
pkin_gen(11) = 0.0;
pkin_gen(13) = 0.0;
