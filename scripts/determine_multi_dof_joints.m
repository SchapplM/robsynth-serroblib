% Bestimme mehrwertige technische Gelenke für die kinematische Ketten
% Suche die geringstmögliche Anzahl dieser Gelenke für alle Modelle in der
% Datenbank
% 
% Siehe auch: MA Zilin He: "Automatic Dimensioning of Drive Trains, Joints
% and Links for Parallel Robots", 2019


% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Initialisierung
roblibpath=fileparts(which('serroblib_path_init.m'));

serroblib_gen_bitarrays(1:7);

%% Durchsuche alle Roboter und prüfe die Kinematikparameter
% Zuordnung der Zahlenwerte in der csv-Tabelle zu den physikalischen Werten
% Die in der mat-Datei abgelegte Binär-Kodierung entspricht der Reihenfolge
for N = 4:6
  fprintf('Prüfe Strukturen mit %d FG\n', N);
  % Alle Roboter aus Datenbank laden
  mdllistfile_Ndof = fullfile(roblibpath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_EEdof0');
  for j = 1:length(l.Names_Ndof)
    RobName = l.Names_Ndof{j};
    fprintf('%d/%d: Prüfe Struktur %s\n', j, length(l.Names_Ndof), RobName);
    RS = serroblib_create_robot_class(RobName);
    RS.gen_testsettings(false, true); % Kinematik-Parameter zufällig setzen
    typestring = RobName(3:3+N-1); % Roboterdaten aus Namen extrahieren
    joint_string = typestring;
    %% Prüfe, welche mehrwertigen Gelenke möglich sind
    % Alle Gelenke durchgehen (siehe Code aus MA Zilin He: optimal_new.m)
    for ii = 1:N
      % Prüfe auf zylindrisches Gelenk
      if ii > 1 && ... % Vorkommen erst ab dem zweiten 1FG-Gelenk auffindbar
          RS.MDH.sigma(ii-1) ~= RS.MDH.sigma(ii) && ... % Dreh- und Schubgelenk nacheinander (Reihenfolge egal)
          RS.MDH.alpha(ii) == 0 ... % Achsen parallel
        % TODO: Ist das so korrekt? Kein Bezug auf Kreuzungsabstand?
        if ~any(joint_string(ii-1) == 'CUS')
          % Nur einsetzen, wenn dadurch nicht vorher schon ein mehrwertiges
          % Gelenk gelöscht wird
          joint_string(ii-1:ii) = '_C'; % Entferne den ersten Eintrag
        end
      end
      % Prüfe auf Kardangelenk
      if ii > 1 && ... % Vorkommen erst ab dem zweiten 1FG-Gelenk auffindbar
          all(RS.MDH.sigma(ii-1:ii) == 0) && ... % Zwei Drehgelenke
          RS.MDH.alpha(ii) == pi/2 && ... % Alle beide senkrecht zueinander
          RS.MDH.a(ii) == 0 && ... % Gelenkachsen schneiden sich in einem Punkt
          RS.MDH.d(ii) == 0
        if ~any(joint_string(ii-1) == 'US')
          if any(joint_string(ii-1) == 'C')
            % Wiederherstellung des einwertigen Gelenks vorher (wird frei
            % durch wegfall des zylindrischen Gelenks)
            joint_string(ii-2) = typestring(ii-2);
          end
          joint_string(ii-1:ii) = '_U'; % Entferne den ersten Eintrag
        end
      end
      % Prüfe auf Kugelgelenk. Da die Kugelgelenk-Prüfung nach dem
      % Kardan-Gelenk kommt, wird es eventuell überschrieben.
      % Das ist so gewollt, da ein Kugelgelenk "besser" ist als ein Kardan-
      % gelenk (als PKM-Beinkette)
      if ii > 2 && ... % Vorkommen erst ab dem dritten 1FG-Gelenk auffindbar
          all(RS.MDH.sigma(ii-2:ii) == 0) && ... % Drei Drehgelenke
          all(RS.MDH.alpha(ii-1:ii) == pi/2) && ... % Alle drei senkrecht zueinander
          all(RS.MDH.a(ii-1:ii) == 0) && ...% Gelenkachsen schneiden sich in einem Punkt
          all(RS.MDH.d(ii-1:ii) == 0)
        if ~any(joint_string(ii-2:ii-1) == 'S')
          % Verhindere, dass ein vorheriges S-Gelenk überschrieben wird.
          % Andere Gelenke wie U/C können überschrieben werden
          if any(joint_string(ii-2) == 'UC')
            % Wiederherstellung des einwertigen Gelenks vorher
            joint_string(ii-3) = typestring(ii-3);
          end
          joint_string(ii-2:ii) = '__S'; % Entferne die ersten beiden Einträge
        end
      end
    end
    joint_string = strrep(joint_string, '_','');
    % Prüfe Wertigkeit der Gelenkfolge
    dof = 3*sum(joint_string=='S') + 2*sum(joint_string=='U') + 2*sum(joint_string=='C') + ...
          1*sum(joint_string=='R') + 1*sum(joint_string=='P');
    if dof ~= N
      error('Die Anzahl der Gelenk-FG stimmt nicht mehr nach der Ersetzung');
    end
    
    if ~strcmp(typestring, joint_string)
      fprintf('Neue Gelenkfolge %s gefunden\n', joint_string);
    else
      % Setze leere Zeichenkette ein, damit man "besondere" Gelenkfolgen
      % schneller erkennt.
      joint_string = '';
    end
    %% Eintragen des Gelenks in die Tabelle
    % Dateinamen der csv-Tabellen bestimmen

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
        csvline_mod{1+8*N+3+9+2} = joint_string;
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
    if ~isempty(joint_string)
      fprintf('\tGelenkfolge auf %s gesetzt.\n', joint_string);
    end
  end
end