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
% EEdof0 [1x9] oder [1x6]
%   Vektor mit beweglichen EE-FG des Roboters (Geschw. und Winkelgeschw. im
%   Basis-KS. Entspricht Vorgabe in der Struktursynthese von Ramirez)
%   Zusätzlich Euler-Winkel Basis-Endeffektor (letzte 3 Komponenten)
%   1="Komponente durch Roboter beeinflussbar"; 0="nicht beeinflussbar"
% 
% Ausgabe:
% mdlname
%   Name des Roboters in der Datenbank
% new
%   true, wenn Roboter neu ist und der Datenbank hinzugefügt wurde.

% Quelle:
% [KhalilKle1986] Khalil, W. and Kleinfinger, J.-F.: "A new geometric
% notation for open and closed-loop robots" (1986) 

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [mdlname, new] = serroblib_add_robot(MDH_struct, EEdof0)

if nargin == 1
  EEdof0 = []; % 6 Leerzeichen
end

if length(EEdof0) == 6
  % Altes Format: Ignoriere die Erweiterung um die Euler-Winkel des EE
  % Setze gewünschte Euler-Winkel auf Winkelgeschwindigkeit (passt bei 1R
  % und 3R; nicht bei 2R). Aber da muss man sowieso die Euler-Winkel nehmen
  EEdof0 = [EEdof0, EEdof0(4:6)];
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
  % Kopfzeile für EE-Trafo erzeugen
  csvline_head1{c} = 'EE-Transformation (phi_N_E)';
  for jj = 1:2
    csvline_head1{c+jj} = '';
  end
  phiNEstr = {'phix_NE', 'phiy_NE', 'phiz_NE'};
  for jj = 1:3
    csvline_head2{c} = phiNEstr{jj}; c=c+1;
  end
  
  % Kopfezeile für EE-FG erzeugen
  csvline_head1{c} = 'EE-FG (Basis-KS)';
  for jj = 1:8
    csvline_head1{c+jj} = '';
  end
  eestr = {'vx0','vy0','vz0','wx0','wy0','wz0','phix0','phiy0','phiz0'};
  for jj = 1:9
    csvline_head2{c} = eestr{jj}; c=c+1;
  end
  % Kopfzeile für weitere Eigenschaften des Roboters
  c = length(csvline_head1)+1;
  csvline_head1{c} = 'Weitere Eigenschaften';
  csvline_head2{c} = 'Positionsbeeinflussendes Gelenk';
  
  % String aus Cell-Array erzeugen
  line_head1 = csvline_head1{1};
  line_head2 = csvline_head2{1};
  for i = 2:length(csvline_head1)
    line_head1 = sprintf('%s;%s', line_head1, csvline_head1{i});
    line_head2 = sprintf('%s;%s', line_head2, csvline_head2{i});
  end
  % Kopfzeile in csv-Tabelle schreiben
  fid = fopen(filepath_csv, 'w');
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
  c = c+1; 
  if strcmp(descr_type{MDH_struct.type(kk)+1}, 'R')
    csvline{c} = sprintf('q%d', kk); %#ok<AGROW>
  else
    csvline{c} = descr_theta{MDH_struct.theta(kk)+1}; %#ok<AGROW>
  end
  c = c+1;
  if strcmp(descr_type{MDH_struct.type(kk)+1}, 'P')
    csvline{c} = sprintf('q%d', kk); %#ok<AGROW>
  else
    csvline{c} = descr_d{MDH_struct.d(kk)+1}; %#ok<AGROW>
  end
  c = c+1; csvline{c} = descr_offset{MDH_struct.offset(kk)+1}; %#ok<AGROW>
end

% Spalten für EE-Transformation
for i = 1:3
  c = c+1; csvline{c} = '?'; %#ok<AGROW>
end

% Spalten für EE-Freiheitsgrade
if ~isempty(EEdof0)
  for i = 1:9
    c = c+1; csvline{c} = sprintf('%d', EEdof0(i)); %#ok<AGROW>
  end
else
  for i = 1:9
    c = c+1; csvline{c} = ''; %#ok<AGROW>
  end
end

% Spalte für Gelenknummer, dass die Position als letztes Beeinflusst
c = c+1; csvline{c} = '?';
%% Zeile für den Roboter finden
% Suche Roboter in den bestehenden csv-Tabellen
[found, idx_direct, num_direct, Name] = serroblib_find_robot(csvline);

% Prüfe, ob der Roboter identisch in der Datenbank ist (ohne Betrachtung
% von Basis-Isomorphismen)
found_ident = serroblib_find_robot(csvline, false);

[found_var, idx_var, ~, Name_var] = serroblib_find_robot(csvline, false, true);

if ~found
  if found_var
    % Hauptmodell des Roboters existiert, aber noch nicht diese Variante.
    % Erhöhe nur die Variantennummer um eins
    index = idx_var(3);
    index_var = idx_var(4) + 1;
  else
    % Roboter existiert noch nicht. Erhöhe die letzte Nummer um eins
    index = num_direct(2)+1;
  end
else
  index = idx_direct(3);
end
% Namen des Robotermodells zusammenstellen.
% Benennungsschema: 
% * S ("Seriell"), 
% * N ("Anzahl FG"), 
% * RRP % ("Gelenkreihenfolge"), 
% * index (Laufende Nummer dieser Konfiguration in allen RRP-Robotern
mdlname = sprintf('S%d%s%d', N, typestring, index);

if found
  if found_ident
    fprintf('serroblib_add_robot: Roboter %s ist schon identisch in der Datenbank.\n', Name);
  else
    fprintf('serroblib_add_robot: Roboter %s ist als Basis-Isomorphismus in der Datenbank.\n', Name);
  end
  new = false;
  new_var = false;
elseif found_var
  new = true;
  new_var = true;
  mdlname_var = sprintf('S%d%s%dV%d', N, typestring, index, index_var);
  fprintf('serroblib_add_robot: Roboter %s ist eine noch unbekannte Variante von %s.\n', mdlname_var, Name_var);
else
  fprintf('serroblib_add_robot: Roboter %s ist noch nicht in der Datenbank.\n', mdlname);
  new = true;
  new_var = false;
end

%% Zeile eines neuen Hauptmodells in Datei hinzufügen
% Zeile ans Ende anhängen
if new && ~new_var
  csvline{1} = mdlname;
  % Cell-Array in csv-Zeile umwandeln (da das Schreiben von Strings nur mit
  % xlswrite funktioniert, dass nicht plattformunabhängig zur Verfügung
  % steht.
  line_robot = mdlname;
  for i = 2:length(csvline)
    line_robot = sprintf('%s;%s', line_robot, csvline{i});
  end

  fid = fopen(filepath_csv, 'a');
  fwrite(fid, [line_robot, newline]);
  fclose(fid);
end

%% Zeile einer neuen Variante in Datei hinzufügen
if new_var
  % Trage Variante in Tabelle an letzter Stelle für das Hauptmodell ein
  if index_var == 1
    mdlname_previous = mdlname;
  else
    mdlname_previous = sprintf('S%d%s%dV%d', N, typestring, index, index_var-1);
  end
  csvline{1} = mdlname_var;
  filename = sprintf('S%d%s.csv', N, typestring);
  filepath_csv = fullfile(repopath, sprintf('mdl_%ddof', N),filename);
  filepath_csv_copy = [filepath_csv, '.copy']; % Kopie der Tabelle zur Bearbeitung
  fid = fopen(filepath_csv);
  fidc = fopen(filepath_csv_copy, 'w');
  tline = fgetl(fid);
  action = false;
  while ischar(tline)
    % Spaltenweise als Cell-Array
    csvline_file = strsplit(tline, ';', 'CollapseDelimiters', false);
    fwrite(fidc, [tline, newline()]);
    if strcmp(csvline_file{1}, mdlname_previous) % Suche Roboternamen (erste Spalte)
      % Neue Zeile für Variante einfügen
      line_new = csvline{1};
      for i = 2:length(csvline)
        line_new = sprintf('%s;%s', line_new, csvline{i});
      end
      fwrite(fidc, [line_new, newline()]);
      action = true;
    end    
    tline = fgetl(fid); % nächste Zeile
  end
  fclose(fid);
  fclose(fidc);
  if ~action
    error('Obwohl eine neue Zeile geschrieben werden sollte, wurde nichts gemacht. Vorheriges Modell %s wurde nicht gefunden. Fehler in Tabellenstruktur von %s', mdlname_previous, filename);
  else
    % Modifizierte Tabelle zurückkopieren
    copyfile(filepath_csv_copy, filepath_csv);
    % Kopie-Tabelle löschen
    delete(filepath_csv_copy);
  end
end

%% Abschluss
% Neu zur mat-Datenbank hinzufügen (wird zum schnelleren Zugriff gepflegt).
% Muss hier erfolgen, damit ein weiterer Aufruf von dieser Funktion den
% vorherigen Roboter auch finden kann (um Nummern korrekt hochzuzählen)
serroblib_gen_bitarrays(N); % Das aktualisiert immer die gesamte Datenbank
