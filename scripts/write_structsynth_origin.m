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

for idx_case = 1:10
  %% Optionen zur Bearbeitung der Tabellen
  flush_data = false;
  set_undef_to_zero = false;
  flush_EEFG_mask = [1 1 1 1 1 1];
  reslist = '';
  switch idx_case
    case 1
      % Name der Datei mit abgespeicherten Namen der kinematischen Ketten
      reslist='structsynth_ser_2T1R';
      % Nummer der Herkunftsspalte (1="Manuell",2="Roboter",3="3T0R-PKM",4="3T1R-PKM")
      idx_oc = 2;
      % Alle Einträge für Roboter mit bestimmten Eigenschaften zurücksetzen
      flush_data = true; flush_Njoint = 3; flush_EEFG = [1 1 0 0 0 1];
    case 2
      reslist='structsynth_ser_3T0R';
      idx_oc = 2;
      flush_data = true; flush_Njoint = 3; flush_EEFG = [1 1 1 0 0 0];
    case 3
      reslist='structsynth_ser_3T1R';
      idx_oc = 2;
      flush_data = true; flush_Njoint = 4; flush_EEFG = [1 1 1 0 0 1];
    case 4
      reslist='structsynth_ser_3T3R';
      idx_oc = 2;
      flush_data = true; flush_Njoint = 6; flush_EEFG = [1 1 1 1 1 1];
    case 5
      reslist='structsynth_ser_3T2R_fixrot';
      idx_oc = 2;
      flush_data = true; flush_Njoint = 5; flush_EEFG = [1 1 1 1 1 0];
    case 6
      % Die Ergebnisse sind aktuell identisch zu structsynth_ser_3T2R_fixrot
      reslist='structsynth_ser_3T2R_varrot';
      idx_oc = 2;
    case 7
      idx_oc = 3; % Spalte für 3T0R-PKM
      flush_data = true;
      flush_Njoint = [4 5]; 
      flush_EEFG = [1 1 1 1 1 1];
      flush_EEFG_mask = [1 1 1 0 0 0]; % Die Rotations-FG sind egal
      reslist = 'structsynth_pkm_3T0R_4J';
    case 8
      idx_oc = 3; % Spalte für 3T0R-PKM
      reslist = 'structsynth_pkm_3T0R_5J';
    case 9
      idx_oc = 4; % Spalte für 3T1R-PKM
      flush_data = true;
      flush_Njoint = 5; 
      flush_EEFG = [1 1 1 1 1 1];
      flush_EEFG_mask = [1 1 1 0 0 0]; % Die Rotations-FG sind egal
      reslist = 'structsynth_pkm_3T1R_5J';
    case 10
      % Letzter Durchgang: Setze alle Felder, in denen bis jetzt ein "?"
      % steht auf 0. Annahme: Alle möglichen Herkünfte der seriellen Ketten
      % wurden bereits vorher aus den Ergebnisliste generiert
      flush_Njoint = 4:5;
      idx_oc = 1:4;
      set_undef_to_zero = true;
      flush_EEFG = [1 1 1 1 1 1];
      flush_EEFG_mask = [0 0 0 0 0 0]; % alle seriellen Ketten bearbeiten
  end
  fprintf('Teil %d: Feststellung der Ergebnisse aus Liste %s\n', idx_case, reslist);
  %% Alle Einträge für bestimmte Roboter zurücksetzen
  % Standardmäßig ist ein Fragezeichen in der csv-Datei gesetzt. Für zu
  % definierende Roboter wird eine "0" gesetzt und anschließend für einige mit
  % einer "1" überschrieben.
  if flush_data || set_undef_to_zero
    for N = flush_Njoint
      % Alle Roboter aus Datenbank laden
      mdllistfile_Ndof = fullfile(roblibpath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
      l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_EEdof0');
      [IndZ, IndB] = serroblib_filter_robots(N, flush_EEFG, flush_EEFG_mask);
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
            if flush_data
              csvline_mod{1+8*N+3+9+1+idx_oc} = '0'; % Setze auf Null
            else % set_undef_to_zero
              % Setze alle Fragezeichen auf Null
              idx = 1+8*N+3+9+1+idx_oc; % Indizes der Spalten zur Modellherkunft
              csvline_mod(idx) = strrep(csvline_mod(idx), '?', '0');
            end
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
          if flush_data
            fprintf('\tWert für Modellherkunft Spalte %d auf 0 gesetzt.\n', idx_oc);
          else
            fprintf('\tWert für Modellherkunft Spalte [%s] auf 0 gesetzt, falls vorher undefiniert.\n', disp_array(idx_oc,'%d'));
          end
        end
      end
    end
  end
  %% Daten aus gespeicherten Ergebnislisten eintragen
  % Nach der Struktursynthese serieller Roboter werden die Ergebnisse als
  % Liste gespeichert. Es wird angenommen, dass die Listen jeweils
  % vollständig sind.
  if isempty(reslist)
    % Falls keine Liste angegeben ist, wird nur eine Löschung vorgenommen
    % (s.o.)
    continue
  end
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
  Names_existing = unique(Names_existing);
  % Speichere die (reduzierte) Liste der tatsächlich in der Datenbank
  % enthaltenen Roboter ab. Durch nachträgliches Entfernen von Ergebnissen
  % der Synthese kann diese Liste kürzer sein.
  writecell(Names_existing(:), [roblist_existing, '.txt']);
  fprintf('Namen der %d aktualisierten Modelle in %s eingetragen.\n', length(Names_existing), roblist_existing);
end
