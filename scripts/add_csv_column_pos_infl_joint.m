% Füge Spalte für Nummer des Pos.-beeinfl. Gelenks zu CSV-Tabellen hinzu.
% Neue Spalte: idx_posinfljoint
% Bearbeitet die Dateien SxPRPR....csv
%
% Vorgehensweise zur Modifikation der csv-Tabellen:
% * Dieses Skript ausführen
% * Anpassung der Funktionen gen_bitarrays, add_robot, csvline2bits, ... an
%   neue Anzahl von Spalten
% * Belegen der Tabelleneinträge mit Werten. 
%   Siehe dazu auch: determine_jointnumber_influence_ee_position.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-04
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Initialisierung
roblibpath=fileparts(which('serroblib_path_init.m'));

%% Durchsuche alle csv-Tabellen und füge Spalten hinzu
for N = 1:7
  mdldir = fullfile(roblibpath, sprintf('mdl_%ddof', N));
  for d = dir(fullfile(mdldir, '*.csv'))'
    % Initialisierung der Dateien
    csvfilepath = fullfile(mdldir, d.name);
    csvfilepath_copy = [csvfilepath, '.copy.csv']; % Kopie der Tabelle zur Bearbeitung
    fid = fopen(csvfilepath);
    fidc = fopen(csvfilepath_copy, 'w');
    tline = fgetl(fid);
    ln = 0; % Zähler für Zeilennummer
    modified = false; % Marker für erfolgte Modifikation der Datei
    error = false; % Marker für Fehler
    while ischar(tline)
      ln = ln+1;
      % Spaltenweise als Cell-Array
      csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
      tline = fgetl(fid); % nächste Zeile
      % Prüfe csv-Zeile. Erkennung anhand der Spaltenzahl, ob Umwandlung
      % bereits erfolgt ist.
      csv_length_old = 1 + 8*N + 3 + 9;
      csv_length_new = 1 + 8*N + 3 + 9 + 1;
      if length(csvline) == csv_length_new
        if modified
          warning('Vorher wurde schon eine Modifikation vorgenommen. Es gab also gemischte Zeilenlängen.');
          error = true;
        end
        break
      elseif length(csvline) ~= csv_length_old
        warning('Zeilenlänge entspricht weder altem, noch neuem Format');
        error = true;
        break;
      end
      % Ab hier kann davon ausgegangen werden, dass das alte Format vorliegt

      % Spalten generieren
      % Siehe serroblib_add_robot.m (Für Konsistenz)
      if ln == 1 % erste Überschriftszeile
        cols_pij = {'Weitere Eigenschaften'};
      elseif ln == 2 % zweite Überschriftszeile
        cols_pij = {'Positionsbeeinflussendes Gelenk'};
      else % Datenzeile mit Roboter
        % Erstmal leer lassen. Das lässt sich an dieser Stelle nicht
        % bestimmen
        cols_pij = {'?'};
      end
      
      % Neue CSV-Zeile erstellen
      csvline_new = {csvline{1:end}, cols_pij{:}}; %#ok<CCAT>
      line_copy = csvline_new{1};
      for i = 2:length(csvline_new)
        line_copy = sprintf('%s;%s', line_copy, csvline_new{i});
      end
      % Zeile in die Dateikopie schreiben
      fwrite(fidc, [line_copy, newline()]);
      modified = true;
    end
    fclose(fid);
    fclose(fidc);
    if modified && ~error
      % Ursprüngliche Datei mit Kopie überschreiben
      copyfile(csvfilepath_copy, csvfilepath);
      delete(csvfilepath_copy);
      fprintf('Datei %s auf neues Format gebracht\n', csvfilepath);
    elseif error
      fprintf('Fehler bei Datei %s\n', csvfilepath);
    else
      fprintf('Keine Änderung in Datei %s\n', csvfilepath);
      delete(csvfilepath_copy);
    end
  end
end