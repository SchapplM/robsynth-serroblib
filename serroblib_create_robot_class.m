% Instanz der Roboterklasse für gegebenen Roboter initialisieren
% 
% Eingabe:
% Name
%   Name des Robotermodells nach dem Schema "SxRRRyy"
%   (x ist die Anzahl der FG, yy die laufende Nummer des Modells dieser
%   Gelenkreihenfolge)
% RobName
%   Name der Roboterparameter entsprechend der ersten Spalte der Tabelle
%   models.csv des jeweiligen Robotermodells
% 
% Ausgabe:
% RS [SerRob]
%   Instanz der SerRob-Klasse mit den Eigenschaften und Funktionen des zu
%   ladenden Roboters

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RS = serroblib_create_robot_class(Name, RobName)

%% Daten laden
N = str2double(Name(2)); % Annahme: Namensschema SxRRR; mit x="Anzahl Gelenk-FG"
repopath=fileparts(which('serroblib_path_init.m'));
mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_phiNE', 'BitArrays_EEdof0');

% Bit-Array für Namen
BA = l.BitArrays_Ndof(strcmp(l.Names_Ndof,Name),:);
if isempty(BA)
  error('Roboter %s ist nicht bekannt', Name);
end
[~, csvbits] = serroblib_bits2csvline(BA);
% Bit-Array für EE-Eigenschaften
BAE = l.BitArrays_EEdof0(strcmp(l.Names_Ndof,Name),:);
[~, EEFG0] = serroblib_bits2csvline_EE(BAE);
% Bit-Array für Endeffektor-Rotation
BAR = l.BitArrays_phiNE(strcmp(l.Names_Ndof,Name),:);
%% Parameter-Struktur erstellen

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


%% Parameter-Zahlenwerte setzen
% In diesem Abschnitt belegte Variablen mit Standardwerten initialisieren:
descr = '';
% Parameter aus Tabelle holen
if nargin > 1 % Falls Name des Parametrierten Modells gegeben
  % Allgemeine Definitionen
  unitmult_angle = pi/180; % Angaben in Tabelle in Grad. Umrechnung in Radiant
  unitmult_dist = 1/1000;
  found = false; % Marker, ob Parametersatz gefunden wurde
  pardat = fullfile(repopath, sprintf('mdl_%ddof', N), Name, 'models.csv');
  if ~exist(pardat, 'file')
    error('Parameterdatei zu %s existiert nicht', Name);
  end
  
  % Suche nach gefordertem Roboter in Parameterdatei
  fid = fopen(pardat);
  tline = fgetl(fid); % csv-Datei zeilenweise einlesen
  while ischar(tline)
    % Spaltenweise als Cell-Array
    csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
    tline = fgetl(fid); % nächste Zeile
    if isempty(csvline) || strcmp(csvline{1}, '')
      continue
    end
    if strcmp(csvline{1}, RobName)
      found = true; % Erste Spalte der Zeile entspricht dem Roboternamen
      break;
    end
  end
  fclose(fid);
  if ~found
    warning('Parameter für %s nicht gefunden', RobName);
  else
    % Daten aus csv-Zeile extrahieren
    c = 3; % erste Drei Zeilen nicht betrachten (sind allgemeine Felder des Roboters)
    descr = sprintf('%s %s', csvline{2}, csvline{3});
    for kk = 1:N % über alle Gelenk-FG
      % Werte aus csv-Zeile herausholen: Einzelne Spalten durchgehen      
      c=c+1; % Gelenktyp (R/P); muss nicht extrahiert werden
      c=c+1; value_beta   = str2double(csvline{c});
      c=c+1; value_b      = str2double(csvline{c});
      c=c+1; value_alpha  = str2double(csvline{c});
      c=c+1; value_a      = str2double(csvline{c});
      c=c+1; value_theta  = str2double(csvline{c});
      c=c+1; value_d      = str2double(csvline{c});
      c=c+1; value_offset = str2double(csvline{c});
      c=c+1; value_qmin   = str2double(csvline{c});
      c=c+1; value_qmax   = str2double(csvline{c});
      c=c+1; value_vmax   = str2double(csvline{c});
      c=c+1; value_sign   = str2double(csvline{c});
      c=c+1; value_qref = str2double(csvline{c});
      
      if PS.sigma(kk) == 0
        unitmult_q = unitmult_angle; % Angabe in Tabelle in deg, Umrechnung in rad
      else
        unitmult_q = unitmult_dist; % Angaben in Tabelle in mm. Umrechnung in m
      end
      % Werte prüfen und eintragen (wenn die Werte ein mit NaN markierter
      % freier Parameter des Robotermodells sind)
      if isnan(PS.beta(kk)), PS.beta(kk)     = unitmult_angle*value_beta; end
      if isnan(PS.b(kk)), PS.b(kk)           = unitmult_dist*value_b; end
      if isnan(PS.alpha(kk)), PS.alpha(kk)   = unitmult_angle*value_alpha; end
      if isnan(PS.a(kk)), PS.a(kk)           = unitmult_dist*value_a; end
      if isnan(PS.theta(kk)), PS.theta(kk)   = unitmult_angle*value_theta; end
      if isnan(PS.d(kk)), PS.d(kk)           = unitmult_dist*value_d; end
      if isnan(PS.offset(kk)), PS.offset(kk) = unitmult_q*value_offset; end
      
      % Werte belegen, wenn sie in der Tabelle nicht gegeben sind
      % Nehme die typischen Einheiten aus der Tabelle (deg, mm) und rechne
      % sie am Ende um
      if isnan(value_qmin)
        if PS.sigma(kk) == 0
          value_qmin = -180;
        else
          value_qmin = -1000;
        end
      end
      if isnan(value_qmax)
        if PS.sigma(kk) == 0
          value_qmax = 180;
        else
          value_qmax = 1000;
        end
      end
      if isnan(value_vmax)
        if PS.sigma(kk) == 0
          value_vmax = 720;
        else
          value_vmax = 1000;
        end
      end
      % Spalte für Achs-Vorzeichen nach Hersteller-Norm (noch nicht vollständig implementiert)
      if value_sign ~= -1 && value_sign ~= 1
        error('Spalte für Vorzeichen hat Wert mit Betrag ungleich eins');
      end
      
      PS.qmin(kk) = unitmult_q*value_qmin;
      PS.qmax(kk) = unitmult_q*value_qmax;
      PS.vmax(kk) = unitmult_q*value_vmax;
      PS.qref(kk) = unitmult_q*value_qref;
    end
  end
end

%% EE-Transformation
phi_N_E = NaN(3,1); % Variable vorbelegen
descr_phi = {0, pi/2, pi, -pi/2, NaN}; % mögliche Werte, die durch Bits kodiert werden
b = 0; % Bit-Zähler
for kk = 1:3 % 3 Euler-Winkel
  % 3 Bits zur Kodierung des aktuellen Wertes extrahieren
  Bit_phi   = bitand( bitshift( BAR, -b), bin2dec('111')); b = b+3;
  % Wert mit den Bits auswählen
  phi_N_E(kk) = descr_phi{Bit_phi+1};
end

%% Klasse initialisieren

% Zahlenwerte der Parameter festlegen.
PS.pkin = []; % sind hier noch gar nicht bekannt. Würde Auswahl eines bestimmten Roboters erfordern

% Klassen-Instanz erstellen
serroblib_addtopath({Name})
RS = SerRob(PS, Name);

% Klassen-Instanz vorbereiten
RS = RS.fill_fcn_handles(false);

% Setze Euler-Winkel für Transformation zum EE-KS
RS.update_EE(zeros(3,1), phi_N_E, uint8(2));

% EE-FG einsetzen
RS.I_EE = logical(EEFG0([1:3,7:9]));
RS.I_EE_Task = RS.I_EE; % Setze EE-FG für Aufgabe (für IK) auf EE-FG des Roboters

RS.update_pkin();

RS.qlim = [PS.qmin(:), PS.qmax(:)];
RS.qDlim = [-PS.vmax(:), PS.vmax(:)];
RS.qref = PS.qref(:);

RS.descr = descr;

if nargin > 1
  % CAD-Modelle initialisieren, falls vorhanden
  cadinidat = fullfile(repopath, sprintf('mdl_%ddof', N), Name, ...
    sprintf('CAD_%s',RobName), sprintf('%s_init_CAD.m', RobName));
  [p,f]=fileparts(cadinidat);
  if exist(cadinidat, 'file')
    addpath(p);
    eval(sprintf('RS = %s(RS, Name, RobName);', f));
    rmpath(p);
  end
end