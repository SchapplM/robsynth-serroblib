% Instanz der Roboterklasse für gegebenen Roboter initialisieren
% 
% Eingabe:
% Name
%   Name des Roboters nach dem Schema "SxRRR"

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function RS = serroblib_create_robot_class(Name)

% Daten laden
N = str2double(Name(2)); % Annahme: Namensschema SxRRR; mit x="Anzahl Gelenk-FG"
repopath=fileparts(which('serroblib_path_init.m'));
mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof');

% Bit-Array für Namen
BA = l.BitArrays_Ndof(strcmp(l.Names_Ndof,Name),:);
[~, csvbits] = serroblib_bits2csvline(BA);

% Parameter-Struktur erstellen

% Mögliche Zustände für MDH-Parameterstruktur bei Eingabe in die Roboterklasse
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
            'sigma', NaN(N,1), ...
            'pkin', [], ...
            'm', NaN(N+1,1), 'mrSges', NaN(N+1,3), 'Ifges', NaN(N+1,6), ...
            'NJ', N, 'NL', N+1, 'NQJ', N);


for kk = 1:N
  PS.beta(kk)  = descr_beta{  csvbits(3+8*(kk-1)) };
  PS.b(kk)     = descr_b{     csvbits(4+8*(kk-1)) };
  PS.alpha(kk) = descr_alpha{ csvbits(5+8*(kk-1)) };  
  PS.a(kk)     = descr_a{     csvbits(6+8*(kk-1)) };
  PS.theta(kk) = descr_theta{ csvbits(7+8*(kk-1)) };
  PS.d(kk)     = descr_d{     csvbits(8+8*(kk-1)) };
end

% Zahlenwerte der Parameter festlegen (TODO)
PS.pkin = []; % TODO

% Klassen-Instanz erstellen
RS = SerRob(PS, Name);

% Klassen-Instanz vorbereiten
RS = RS.fill_fcn_handles(false);

