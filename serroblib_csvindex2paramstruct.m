% Wandle die CSV-Index-Zeile in eine Parameterstruktur um
% 
% Eingabe:
% N
%   Anzahl der Gelenke
% csvbits
%   Spaltenweise Indizes für die Daten in csvline (Zeile der csv-Tabelle)
%   bezogen auf die möglichen Werte für jeden Eintrag. 
%   Siehe serroblib_bits2csvline
% 
% Ausgabe:
% PS
%   Parameter-Struktur (Eingabe für SerRob-Klasse)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2023-03
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function PS = serroblib_csvindex2paramstruct(N, csvbits)
% Mögliche Zustände für MDH-Parameterstruktur bei Eingabe in die
% Roboterklasse (siehe auch serroblib_csvline2bits.m)
descr_type = {0, 1};
descr_beta = {0, pi/2, pi, -pi/2, NaN};
descr_b = {0, NaN};
descr_alpha = {0, pi/2, pi, -pi/2, NaN};
descr_a = {0, NaN};
descr_theta = {0, pi/2, pi, -pi/2, NaN};
descr_d = {0, NaN};
descr_offset = {0, pi/2, pi, -pi/2, NaN};

% Parameter-Struktur für Eingabe in Roboterklasse
% Siehe auch: serroblib_generate_mapleinput.m (zur Verwendung von csvbits)
PS = struct('beta',  NaN(N,1), 'b', NaN(N,1), ...
            'alpha', NaN(N,1), 'a', NaN(N,1), ...
            'theta', NaN(N,1), 'd', NaN(N,1), ...
            'sigma', NaN(N,1), 'offset', NaN(N,1), ...
            'pkin', [], 'v', uint8(0:N-1)', ...
            'mu', ones(N,1), ...
            'm', NaN(N+1,1), 'mrSges', NaN(N+1,3), 'Ifges', NaN(N+1,6), ...
            'NJ', N, 'NL', N+1, 'NQJ', N, ...
            'qmin', NaN(N,1), 'qmax', NaN(N,1), 'vmax', NaN(N,1), 'qref', NaN(N,1));

for kk = 1:N
  PS.sigma(kk) = descr_type{  csvbits(2+8*(kk-1)) };
  PS.beta(kk)  = descr_beta{  csvbits(3+8*(kk-1)) };
  PS.b(kk)     = descr_b{     csvbits(4+8*(kk-1)) };
  PS.alpha(kk) = descr_alpha{ csvbits(5+8*(kk-1)) };  
  PS.a(kk)     = descr_a{     csvbits(6+8*(kk-1)) };
  PS.theta(kk) = descr_theta{ csvbits(7+8*(kk-1)) };
  PS.d(kk)     = descr_d{     csvbits(8+8*(kk-1)) };
  PS.offset(kk)= descr_offset{csvbits(9+8*(kk-1)) };
end
