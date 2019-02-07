% Füge Spalten für erweiterte Darstellung der EE-FG zu CSV-Tabellen hinzu.
% Neue Spalten: phi_N_E (Platzhalter konstante Trafo), phi_0_E
% (Orientierungs-FG in XYZ-Euler-Winkeln)
% Bearbeitet die Dateien SxPRPR....csv
% 
% Siehe auch: correct_phi_N_E.m (zum Einsetzen von Werten in die Spalten
% für phi_N_E)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
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
      csv_length_old = 1 + 8*N + 6;
      csv_length_new = 1 + 8*N + 3 + 6 + 3;
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
      % cols_phiNE: Rotation vom letzten Körper-KS (N) zum EE-KS (E)
      % cols_phi0E: Rotation von Basis zum Endeffektor (als Euler-XYZ)
      % Siehe serroblib_add_robot.m (Für Konsistenz)
      if ln == 1 % erste Überschriftszeile
        cols_phiNE = {'EE-Transformation (phi_N_E)', '', ''};
        cols_phi0E = {'', '', ''};
      elseif ln == 2 % zweite Überschriftszeile
        cols_phiNE = {'phix_NE', 'phiy_NE', 'phiz_NE'};
        cols_phi0E = {'phix0', 'phiy0', 'phiz0'};
      else % Datenzeile mit Roboter
        % Erstmal leer lassen. Das lässt sich an dieser Stelle nicht
        % bestimmen
        cols_phiNE = {'?', '?', '?'};
        % Rotation zum EE aus Winkelgeschwindigkeit zum EE bestimmen
        [~, BAE] = serroblib_csvline2bits(csvline); % Hier wurde noch die alte Version benutzt, für die alte csv-Version
        if     all(bitget(BAE,4:6) == [0 0 1])
          cols_phi0E = {'0', '0', '1'};
        elseif all(bitget(BAE,4:6) == [1 1 1])
          cols_phi0E = {'1', '1', '1'};
        elseif all(bitget(BAE,4:6) == [0 0 0])
          cols_phi0E = {'0', '0', '0'};
        else
          warning('Nicht bestimmbar');
          error = true;
          break;
        end
      end
      
      % Neue CSV-Zeile erstellen
      csvline_new = { ...
        csvline{(1):(1 + 8*N)}, ...
        cols_phiNE{:}, ...
        csvline{(1 + 8*N + 1):end}, ...
        cols_phi0E{:}}; %#ok<CCAT>
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