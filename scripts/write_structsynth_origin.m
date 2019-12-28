% Schreibe die Modellherkunft der kinematischen Ketten in die Datenbank
% Dadurch kann im Nachhinein festgestellt werden, wie eine serielle Kette
% in die Datenbank gelangt ist

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Initialisierung
roblibpath=fileparts(which('serroblib_path_init.m'));
robot_list_dir = fullfile(roblibpath, 'synthesis_result_lists');
serroblib_gen_bitarrays(1:7);

for idx_case = 1:4
  %% Optionen zur Bearbeitung der Tabellen
  flush_data = false;
  switch idx_case
    case 1
      % Name der Datei mit abgespeicherten Namen der kinematischen Ketten
      reslist='structsynth_ser_2T1R';
      % Nummer der Herkunftsspalte (1="Roboter",2="3T0R-PKM",3="3T1R-PKM")
      idx_oc = 1;
      % Alle Einträge für Roboter mit bestimmten Eigenschaften zurücksetzen
      flush_data = true; flush_Njoint = 3; flush_EEFG = [1 1 0 0 0 1];
    case 2
      reslist='structsynth_ser_3T0R';
      idx_oc = 1;
      flush_data = true; flush_Njoint = 3; flush_EEFG = [1 1 1 0 0 0];
    case 3
      reslist='structsynth_ser_3T1R';
      idx_oc = 1;
      flush_data = true; flush_Njoint = 4; flush_EEFG = [1 1 1 0 0 1];
    case 4
      reslist='structsynth_ser_3T3R';
      idx_oc = 1;
      flush_data = true; flush_Njoint = 6; flush_EEFG = [1 1 1 1 1 1];
  end
  fprintf('Teil %d: Feststellung der Ergebnisse aus Liste %s\n', idx_case, reslist);
  %% Alle Einträge für bestimmte Roboter zurücksetzen
  % Standardmäßig ist ein Fragezeichen in der csv-Datei gesetzt. Für zu
  % definierende Roboter wird eine "0" gesetzt und anschließend für einige mit
  % einer "1" überschrieben.
  if flush_data
    EE_FG = flush_EEFG;
    for N = flush_Njoint
      % Alle Roboter aus Datenbank laden
      mdllistfile_Ndof = fullfile(roblibpath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
      l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_EEdof0');
      [IndZ, IndB] = serroblib_filter_robots(N, EE_FG, EE_FG);
      k = 0;
      for j = IndZ'
        k = k+1;
        RobName = l.Names_Ndof{j};
        fprintf('%d/%d: Prüfe Struktur %s\n', k, length(IndZ), RobName);
        % Aktualisieren in Tabelle
        typestring = RobName(3:3+N-1); % Roboterdaten aus Namen extrahieren
        filename = sprintf('S%d%s.csv', N, typestring);
        filepath_csv = fullfile(roblibpath, sprintf('mdl_%ddof', N),filename);
        filepath_csv_copy = [filepath_csv, '.copy']; % Kopie der Tabelle zur Bearbeitung
        % CSV-Tabelle zeilenweise durchgehen
        fid = fopen(filepath_csv);
        fidc = fopen(filepath_csv_copy, 'w');
        tline = fgetl(fid);
        found = false;
        while ischar(tline)
          % Spaltenweise als Cell-Array
          csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
          if strcmp(csvline{1}, RobName) % Suche Roboternamen (erste Spalte)
            % Zu modifizierenden Roboter gefunden. Diese Zeile verändert in Dateikopie schreiben
            found = true;
            % Zeile modifizieren
            csvline_mod = csvline;
            csvline_mod{1+8*N+3+9+1+idx_oc} = '0'; % Setze auf Null
            % modifizierte csv-Zeile in Textzeile umwandeln
            line_copy = csvline_mod{1};
            for i = 2:length(csvline_mod)
              line_copy = sprintf('%s;%s', line_copy, csvline_mod{i});
            end
            fwrite(fidc, [line_copy, newline()]);
          else
            % ursprüngliche Zeile in die Dateikopie schreiben: Überschriften,
            % andere Roboter
            fwrite(fidc, [tline, newline()]);
          end
          tline = fgetl(fid); % nächste Zeile
        end
        fclose(fid);
        fclose(fidc);
        if ~found
          % Womöglich wurde das Modell schon gelöscht, weil das Ergebnis
          % der Struktursynthese schlecht war.
          warning('Zu bearbeitendes Modell %s nicht in %s gefunden', ...
            RobName, filepath_csv);
        else
          % Modifizierte Tabelle zurückkopieren
          copyfile(filepath_csv_copy, filepath_csv);
        end
        % Kopie-Tabelle löschen
        delete(filepath_csv_copy);
        if found
          fprintf('\tWert für Modellherkunft Spalte %d auf 0 gesetzt.\n', idx_oc);
        end
      end
    end
  end
  %% Daten aus gespeicherten Ergebnislisten eintragen
  % Nach der Struktursynthese serieller Roboter werden die Ergebnisse als
  % Liste gespeichert. Es wird angenommen, dass die Listen jeweils
  % vollständig sind.
  roblist = fullfile(robot_list_dir, reslist);
  roblist_existing = fullfile(robot_list_dir, [reslist, '_in_DB']);
  Names = readcell([roblist, '.txt']);
  Names_existing = {}; % reduzierte Liste von tatsächlich existierenden Robotern in der Datenbank
  for j = 1:length(Names) % Robternamen aus Ergebnisliste durchgehen
    RobName = Names{j};
    fprintf('%d/%d: Aktualisiere Herkunft der Struktur %s\n', j, length(Names), RobName);
    % Aktualisieren in Tabelle
    N = str2double(RobName(2));
    typestring = RobName(3:3+N-1); % Roboterdaten aus Namen extrahieren
    filename = sprintf('S%d%s.csv', N, typestring);
    filepath_csv = fullfile(roblibpath, sprintf('mdl_%ddof', N),filename);
    filepath_csv_copy = [filepath_csv, '.copy']; % Kopie der Tabelle zur Bearbeitung
    % CSV-Tabelle zeilenweise durchgehen
    fid = fopen(filepath_csv);
    fidc = fopen(filepath_csv_copy, 'w');
    tline = fgetl(fid);
    found = false;
    while ischar(tline)
      % Spaltenweise als Cell-Array
      csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
      if strcmp(csvline{1}, RobName) % Suche Roboternamen (erste Spalte)
        % Zu modifizierenden Roboter gefunden. Diese Zeile verändert in Dateikopie schreiben
        found = true;
        % Zeile modifizieren
        csvline_mod = csvline;
        csvline_mod{1+8*N+3+9+1+idx_oc} = '1'; % Setze auf Eins
        % modifizierte csv-Zeile in Textzeile umwandeln
        line_copy = csvline_mod{1};
        for i = 2:length(csvline_mod)
          line_copy = sprintf('%s;%s', line_copy, csvline_mod{i});
        end
        fwrite(fidc, [line_copy, newline()]);
      else
        % ursprüngliche Zeile in die Dateikopie schreiben: Überschriften,
        % andere Roboter
        fwrite(fidc, [tline, newline()]);
      end
      tline = fgetl(fid); % nächste Zeile
    end
    fclose(fid);
    fclose(fidc);
    if ~found
      warning('Zu bearbeitendes Modell %s nicht in %s gefunden', ...
        RobName, filepath_csv);
    else
      % Modifizierte Tabelle zurückkopieren
      copyfile(filepath_csv_copy, filepath_csv);
    end
    % Kopie-Tabelle löschen
    delete(filepath_csv_copy);
    if found
      Names_existing = [Names_existing(:)', {RobName}];
      fprintf('\tWert für Modellherkunft Spalte %d auf 1 gesetzt.\n', idx_oc);
    end
  end
  % Speichere die (reduzierte) Liste der tatsächlich in der Datenbank
  % enthaltenen Roboter ab. Durch nachträgliches Entfernen von Ergebnissen
  % der Synthese kann diese Liste kürzer sein.
  writecell(unique(Names_existing(:)), [roblist_existing, '.txt']);
end