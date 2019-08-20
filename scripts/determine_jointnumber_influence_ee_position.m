% Bestimme die Gelenk-Nummer, die die Position des Endeffektor-Segmentes
% beeinflusst.
% Dadurch kann entschieden werden, ob die serielle Kette als PKM-Beinkette
% geeignet ist und welche Modellierung dort benutzt werden muss.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-04
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Initialisierung
roblibpath=fileparts(which('serroblib_path_init.m'));

serroblib_gen_bitarrays(1:7);

%% Durchsuche alle Roboter und prüfe die Kinematikparameter
% Zuordnung der Zahlenwerte in der csv-Tabelle zu den physikalischen Werten
% Die in der mat-Datei abgelegte Binär-Kodierung entspricht der Reihenfolge
for N = 2:7
  fprintf('Prüfe Strukturen mit %d FG\n', N);
  % Alle Roboter aus Datenbank laden
  mdllistfile_Ndof = fullfile(roblibpath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_EEdof0');
  for j = 1:length(l.Names_Ndof)
    RobName = l.Names_Ndof{j};
    fprintf('%d/%d: Prüfe Struktur %s\n', j, length(l.Names_Ndof), RobName);
    RS = serroblib_create_robot_class(RobName);
    RS.gen_testsettings(false, true); % Kinematik-Parameter zufällig setzen
    %% Prüfe, welches Gelenk die EE-Position beeinflusst
    % Prüfe, für jedes Gelenk beginnend von hinten, ob es die EE-Position
    % beeinflusst
    idx_joint_posinfl = 0;
    q = rand(N,1);
    J = RS.jacobig(q);
    for i_joint = N:-1:1
      % Wenn das Gelenk die Position beeinflussen kann, ist die
      % entsprechende Jacobi-Spalte ungleich Null
      if any(J(1:3,i_joint))
        % Dieses Gelenk ist das erste (vom EE gezählt), dass die
        % EE-Position beeinflusst.
        idx_joint_posinfl = i_joint;
      	break;
      end
    end
    %% Eintragen des Gelenks in die Tabelle
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
    while ischar(tline)
      % Spaltenweise als Cell-Array
      csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
      if strcmp(csvline{1}, RobName) % Suche Roboternamen (erste Spalte)
        % Zu modifizierenden Roboter gefunden. Diese Zeile verändert in Dateikopie schreiben
        found = true;
        % Zeile modifizieren
        csvline_mod = csvline;
        csvline_mod{1+8*N+3+9+1} = sprintf('%d', idx_joint_posinfl);
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
    % Modifizierte Tabelle zurückkopieren
    copyfile(filepath_csv_copy, filepath_csv);
    % Kopie-Tabelle löschen
    delete(filepath_csv_copy);
    fprintf('\tWert Pos.-beeinfl. Gelenk = %d gesetzt.\n', idx_joint_posinfl);
  end
end