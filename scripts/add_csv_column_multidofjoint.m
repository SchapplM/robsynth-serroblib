% Füge Spalte für mehrwertige Gelenke zu CSV-Tabellen hinzu.
% In der Spalte ist eingetragen, ob eine Kugel- oder Kardangelenk möglich
% ist (bzw. die genaue Gelenkfolge wie z.B "UPS")
%
% Vorgehensweise zur Modifikation der csv-Tabellen:
% * Dieses Skript ausführen
% * Anpassung der Funktionen gen_bitarrays, add_robot, csvline2bits, ... an
%   neue Anzahl von Spalten
% * Belegen der Tabelleneinträge mit Werten.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Initialisierung
roblibpath=fileparts(which('serroblib_path_init.m'));

%% Durchsuche alle csv-Tabellen und füge Spalten hinzu
for N = 2:7
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
      csv_length_old = 1 + 8*N + 3 + 9 + 1 + 4;
      csv_length_new = csv_length_old + 2; % 2 neue Spalten
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

      % Spalten generieren. Füge letzte drei Spalten der bisherigen Zeile
      % rechts von der hinzuzufügenden Spalte an.
      % Siehe serroblib_add_robot.m (Für Konsistenz)
      if ln == 1 % erste Überschriftszeile
        cols_mj1 = {''}; % Unterhalb der Überschrift "Weitere Eigenschaften"
        cols_mj2 = {''}; % Unterhalb der Überschrift "Herkunft Struktursynthese"
      elseif ln == 2 % zweite Überschriftszeile
        cols_mj1 = {'Gelenkfolge'};
        cols_mj2 = {'Gelenkfolge'};
      else % Datenzeile mit Roboter
        % Erstmal auf Null setzen. Muss manuell eingetragen werden.
        cols_mj1 = {''}; % Hier stehen später Einträge wie "UPS"
        cols_mj2 = {'?'}; % Hier steht eine 1, wenn Gelenkfolge aus automatischer Generierung
      end
      
      % Neue CSV-Zeile erstellen (am Ende von Rubrik "weitere
      % Eigenschaften" und "Erkunft" jeweils eintragen)
      csvline_new = {csvline{1:end-4}, cols_mj1{:}, csvline{end-3:end}, cols_mj2{:}}; %#ok<CCAT>
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