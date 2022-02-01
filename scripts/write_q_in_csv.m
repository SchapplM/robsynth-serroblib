% Schreibe die Werte für q in die CSV-Tabellen in die Spalte für d/theta. 
% Dient zur Korrektur bei versehentlichem Entfernen in anderem Skript. 
% Hat keinen Einfluss auf die Funktionalität. Zur besseren Übersicht.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Initialisierung
roblibpath=fileparts(which('serroblib_path_init.m'));
serroblib_gen_bitarrays(1:7);

%% Durchsuche alle Roboter und prüfe die Einträge
for N = 1:7
  fprintf('Prüfe Strukturen mit %d FG\n', N);
  % Alle Roboter aus Datenbank laden
  mdllistfile_Ndof = fullfile(roblibpath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_EEdof0', 'AdditionalInfo');
  for j = 1:length(l.Names_Ndof)
    RobName = l.Names_Ndof{j};
    fprintf('%d/%d: Prüfe Struktur %s\n', j, length(l.Names_Ndof), RobName);

    %% Eintragen der EE-FG in die Tabelle
    % Dateinamen der csv-Tabellen bestimmen
    typestring = RobName(3:3+N-1); % Roboterdaten aus Namen extrahieren
    filename = sprintf('S%d%s.csv', N, typestring);
    filepath_csv = fullfile(roblibpath, sprintf('mdl_%ddof', N),filename);
    filepath_csv_copy = [filepath_csv, '.copy']; % Kopie der Tabelle zur Bearbeitung
    % CSV-Tabelle zeilenweise durchgehen
    fid = fopen(filepath_csv);
    fidc = fopen(filepath_csv_copy, 'w');
    tline = fgetl(fid);
    found = false;
    changed = false;
    while ischar(tline)
      % Spaltenweise als Cell-Array
      csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
      if strcmp(csvline{1}, RobName) % Suche Roboternamen (erste Spalte)
        % Zu modifizierenden Roboter gefunden. Diese Zeile verändert in Dateikopie schreiben
        found = true;
        % Zeile modifizieren: Eintrag für q1, q2, ... schreiben
        csvline_mod = csvline;
        for jj = 1:N
          if typestring(jj) == 'P'
            idxoffset = 7; % Zur Indizierung von Spalte d
          else
            idxoffset = 6; % Zur Indizierung von Spalte theta
          end
          csvline_mod{1+8*(jj-1)+idxoffset} = sprintf('q%d', jj);
        end
        % Prüfe, ob modifizierte Zeile anders ist
        I_change = ~strcmp(csvline_mod, csvline);
        if any(I_change)
          fprintf('\tZeile %s wurde geändert. Spalten [%s]: {%s} -> {%s}\n', ...
            RobName, disp_array(find(I_change), '%d'), disp_array(...
            csvline(I_change),'%s'), disp_array(csvline_mod(I_change), '%s'));
          changed = true;
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
      warning('Zu bearbeitendes Modell %s nicht in %s gefunden', ...
        RobName, filepath_csv);
      return
    end
    if changed % Es wurden Daten geändert. Nur in diesem Fall DB bearbeiten
      % Modifizierte Tabelle zurückkopieren
      copyfile(filepath_csv_copy, filepath_csv);
    end
    % Kopie-Tabelle löschen
    delete(filepath_csv_copy);
    if changed
      fprintf('\tEintrag in Datenbank geändert.\n');
    end
  end
end
