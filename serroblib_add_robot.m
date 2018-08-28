% Füge ein Robotermodell zur Datenbank (csv-Tabelle) hinzu
% Das Robotermodell wird in die csv-Tabelle für den Gelenktyp (z.B. RPR)
% geschrieben
% 
% Eingabe:
% MDH_struct
%   Struktur mit Feldern für MDH-Parameter
%   Die Felder enthalten nicht die Information selbst (z.B. "Pi/2", sondern
%   einen Index (beginnend bei 0), der auf die Information zeigt)
%   Siehe dazu auch serroblib_bits2csvline.m, serroblib_csvline2bits.m
%   type   Dreh- oder Schub-Gelenk
%   beta   Rotation z
%   b      Translation z
%   alpha  Rotation x
%   a      Translation x
%   theta  Rotation z (Gelenkkoordinate)
%   d      Translation z (Gelenkkoordinate)
%   offset Gelenk-Offset
% EEdof0
%   Vektor mit beweglichen EE-FG des Roboters (Geschw. und Winkelgeschw. im
%   Basis-KS. Entspricht Vorgabe in der Struktursynthese von Ramirez)
%   1="Komponente durch Roboter beeinflussbar"; 0="nicht beeinflussbar"
% 
% Ausgabe:
% mdlname
%   Name des Roboters in der Datenbank

% Quelle:
% [KhalilKle1986] Khalil, W. and Kleinfinger, J.-F.: "A new geometric
% notation for open and closed-loop robots" (1986) 

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function mdlname = serroblib_add_robot(MDH_struct, EEdof0)

if nargin == 1
  EEdof0 = []; % 6 Leerzeichen
end

% Anzahl Gelenke
N = length(MDH_struct.type);

% Typkennzeichnung des Roboters erzeugen (Gelenkreihenfolge)
typestring = char(zeros(1,N));
for i = 1:N
  if MDH_struct.type(i) == 0
    typestring(i) = 'R';
  else
    typestring(i) = 'P';
  end
end

%% Modell in csv-Tabelle eintragen
% Die Tabelle enthält alle möglichen MDH-Parameter für diesen Robotertyp
% (z.B. alle Varianten für RRR)
filename = sprintf('S%d%s.csv', N, typestring);
repopath=fileparts(which('serroblib_path_init.m'));
filepath_csv = fullfile(repopath, sprintf('mdl_%ddof', N), filename);

%% Kopfzeile für Tabelle zusammenstellen
if ~exist(filepath_csv, 'file')
  % Datei existiert nicht. Erstelle Kopfzeile
  csvline_head1 = {'Name'};
  csvline_head2 = {''};
  c = 2; % Spalten-Zähler
  % Kopfzeile für alle Gelenke erstellen
  for kk = 1:N
    csvline_head1{c} = sprintf('Gelenk %d', kk);
    for jj = 1:7
      csvline_head1{c+jj} = '';
    end
    csvline_head2{c} = 'Typ'; c = c+1;
    csvline_head2{c} = 'beta'; c = c+1;
    csvline_head2{c} = 'b'; c = c+1;
    csvline_head2{c} = 'alpha'; c = c+1;
    csvline_head2{c} = 'a'; c = c+1;
    csvline_head2{c} = 'theta'; c = c+1;
    csvline_head2{c} = 'd'; c = c+1;
    csvline_head2{c} = 'offset'; c = c+1;
  end
  % Kopfezeile für EE-FG erzeugen
  csvline_head1{c} = 'EE-FG (Basis-KS)';
  for jj = 1:5
    csvline_head1{c+jj} = '';
  end
  eestr = {'vx0','vy0','vz0','wx0','wy0','wz0'};
  for jj = 1:6
    csvline_head2{c} = eestr{jj}; c=c+1;
  end
  % String aus Cell-Array erzeugen
  line_head1 = csvline_head1{1};
  line_head2 = csvline_head2{1};
  for i = 2:length(csvline_head1)
    line_head1 = sprintf('%s,%s', line_head1, csvline_head1{i});
    line_head2 = sprintf('%s,%s', line_head2, csvline_head2{i});
  end
  % Kopfzeile in csv-Tabelle schreiben
  fid = fopen(filepath_csv, 'a');
  fwrite(fid, [line_head1, newline]);
  fwrite(fid, [line_head2, newline]);
  fclose(fid);
end
%% csv-Tabellenzeile zusammenstellen
% Zeile für die Roboterkinematik in der Tabelle
csvline = {};
csvline{1} = 'Name unknown'; % wird später belegt
c = 1;

for kk = 1:N
  % Liste von möglichen Zuständen für die einzelnen Parameter. Die
  % Reihenfolge dieser Zustände entspricht den Indizes in anderen
  % Funktionen
  descr_type = {'R', 'P'};
  descr_beta = {'0', 'pi/2', 'pi', '-pi/2', sprintf('beta%d',kk)};
  descr_b = {'0', sprintf('b%d',kk)};
  descr_alpha = {'0', 'pi/2', 'pi', '-pi/2', sprintf('alpha%d',kk)};
  descr_a = {'0', sprintf('a%d',kk)};
  descr_theta = {'0', 'pi/2', 'pi', '-pi/2', sprintf('theta%d',kk)};
  descr_d = {'0', sprintf('d%d',kk)};
  descr_offset = {'0', 'pi/2', 'pi', '-pi/2', sprintf('offset%d',kk)};
  
  % CSV-Spalten für aktuelles Gelenk "kk" zusammenstellen
  % (Indizes auf Daten kommen aus MDH_struct)
  c = c+1; csvline{c} = descr_type{MDH_struct.type(kk)+1}; %#ok<AGROW>
  c = c+1; csvline{c} = descr_beta{MDH_struct.beta(kk)+1}; %#ok<AGROW>
  c = c+1; csvline{c} = descr_b{MDH_struct.b(kk)+1}; %#ok<AGROW>
  c = c+1; csvline{c} = descr_alpha{MDH_struct.alpha(kk)+1}; %#ok<AGROW>
  c = c+1; csvline{c} = descr_a{MDH_struct.a(kk)+1}; %#ok<AGROW>
  c = c+1; csvline{c} = descr_theta{MDH_struct.theta(kk)+1}; %#ok<AGROW>
  c = c+1; csvline{c} = descr_d{MDH_struct.d(kk)+1}; %#ok<AGROW>
  c = c+1; csvline{c} = descr_offset{MDH_struct.offset(kk)+1}; %#ok<AGROW>
end

% Spalten für EE-Freiheitsgrade
if ~isempty(EEdof0)
  for i = 1:6
    c = c+1; csvline{c} = EEdof0(i); %#ok<AGROW>
  end
else
  for i = 1:6
    c = c+1; csvline{c} = ''; %#ok<AGROW>
  end
end

%% Zeile für den Roboter finden
% Suche Roboter in den bestehenden csv-Tabellen
[found, index, num] = serroblib_find_robot(csvline);
if ~found
  % Roboter existiert noch nicht. Erhöhe die letzte Nummer um eins
  index = num+1;
end
% Namen des Robotermodells zusammenstellen.
% Benennungsschema: 
% * S ("Seriell"), 
% * N ("Anzahl FG"), 
% * RRP % ("Gelenkreihenfolge"), 
% * index (Laufende Nummer dieser Konfiguration in allen RRP-Robotern
mdlname = sprintf('S%d%s%d', N, typestring, index);

csvline{1} = mdlname;
if found
  fprintf('serroblib_add_robot: Roboter %s ist schon in der Datenbank.\n', mdlname);
  return;
else
  fprintf('serroblib_add_robot: Roboter %s ist noch nicht in der Datenbank.\n', mdlname);
end

%% Zeile in Datei hinzufügen

% Cell-Array in csv-Zeile umwandeln (da das Schreiben von Strings nur mit
% xlswrite funktioniert, dass nicht plattformunabhängig zur Verfügung
% steht.
line_robot = mdlname;
for i = 2:length(csvline)
  line_robot = sprintf('%s,%s', line_robot, csvline{i});
end

fid = fopen(filepath_csv, 'a');
fwrite(fid, [line_robot, newline]);
fclose(fid);

% Neu zur mat-Datenbank hinzufügen (wird zum schnelleren Zugriff gepflegt).
% Muss hier erfolgen, damit ein weiterer Aufruf von dieser Funktion den
% vorherigen Roboter auch finden kann (um Nummern korrekt hochzuzählen)
serroblib_gen_bitarrays(N); % Das aktualisiert immer die gesamte Datenbank

