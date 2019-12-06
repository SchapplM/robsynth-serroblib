% Stelle die EE-FG der Strukturen in der Datenbank fest und schreibe diese
% in die entsprechenden Tabellen.
% Diese Eigenschaft wird normalerweise bei der Struktursynthese erzeugt.
% Bei händisch erstellten Strukturen fehlt sie aber eventuell.
% Mit diesem Skript kann die Eigenschaft auch geprüft werden.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Initialisierung
roblibpath=fileparts(which('serroblib_path_init.m'));
serroblib_gen_bitarrays(1:7);

% Einstellungen
usr_dryrun = false; % Nur Anzeigen, was gemacht werden würde, nichts schreiben
usr_overwrite = true; % Überschreibe auch bestehende Einträge. Sinnvoll, wenn die Spalten vorher leer sind
usr_abortonerror = true; % Bei irgendeinem Fehler anhalten
only_look_at_robot = {}; % Nur eine Liste namentlich genannter Roboter bearbeiten; z.B. S5PRPRR4 

%% Durchsuche alle Roboter und prüfe die Kinematikparameter
% Zuordnung der Zahlenwerte in der csv-Tabelle zu den physikalischen Werten
% Die in der mat-Datei abgelegte Binär-Kodierung entspricht der Reihenfolge
for N = 1:7
  fprintf('Prüfe Strukturen mit %d FG\n', N);
  % Alle Roboter aus Datenbank laden
  mdllistfile_Ndof = fullfile(roblibpath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_EEdof0');
  for j = 1:length(l.Names_Ndof)
    RobName = l.Names_Ndof{j};
    if ~isempty(only_look_at_robot) && ~any(strcmp(only_look_at_robot, RobName))
      continue
    end
    fprintf('%d/%d: Prüfe Struktur %s\n', j, length(l.Names_Ndof), RobName);
    RS = serroblib_create_robot_class(RobName);
    try
     RS.gen_testsettings(false, true);
    catch
      warning('Fehler für Modell %s', RobName);
      if usr_abortonerror
        return
      else
        continue
      end
    end
    if any(isnan(RS.phi_N_E))
      warning('EE-Transformation N-E wurde noch nicht gesetzt. Setze vorerst zu Null');
      RS.update_EE(zeros(3,1), zeros(3,1));
    end
    
    %% Prüfe, welche EE-FG die Struktur hat
    % Prüfe, für jedes Gelenk beginnend von hinten, ob es die EE-Position
    % beeinflusst
    q = rand(N,1);
    qD = rand(N,1);
    Jg = RS.jacobig(q);
    v = Jg*qD;
    Ja = RS.jacobia(q);
    RS.phi_N_E;
    RS.T_N_E;
    T_E = RS.fkineEE(q);
    xD = Ja*qD;
    if any(isnan(xD))
      warning('%s: Euler-Geschwindigkeit kann wegen Singularität nicht berechnet werden', RobName);
      if usr_abortonerror
        return
      else
        continue
      end
    end
    % Bestimme die EE-FG aufgrund der Existenz der
    % Geschwindigkeits-Komponenten
    I_EE_0 = abs([v;xD(4:6)]) > 1e-9;
    
    % Sonderfall 3T2R: Setze die letzte Rotation (nur Euler-Winkel) auf
    % Null (Konvention zum Abspeichern und kennzeichnen)
    if rank(Ja(4:6,:)) == 2 && rank(Ja(4:5,:)) == 2
      I_EE_0(9) = 0;
    end
    
    csvline_EE_FG_calc = cell(1,9);
    csvline_EE_FG_calc(I_EE_0) = {'1'};
    csvline_EE_FG_calc(~I_EE_0) = {'0'};
    
    if serroblib_csvline2bits_EE(csvline_EE_FG_calc) ~= l.BitArrays_EEdof0(j)
      warning('Der Roboter hat in der DB andere EE-FG (%s) als aus der Berechnung (%s)', ...
        dec2bin(l.BitArrays_EEdof0(j),9), dec2bin(serroblib_csvline2bits_EE(csvline_EE_FG_calc),9));
    elseif ~usr_overwrite
      fprintf('Die EE-FG von %s sind bereits richtig in DB (%s). Nichts unternehmen\n', ...
        RobName, dec2bin(l.BitArrays_EEdof0(j),9));
      continue
    else
      fprintf('Die EE-FG von %s sind bereits richtig in DB (%s). Überschreibe trotzdem\n', ...
        RobName, dec2bin(l.BitArrays_EEdof0(j),9));
    end
    if usr_dryrun
      continue
    end
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
    while ischar(tline)
      % Spaltenweise als Cell-Array
      csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
      if strcmp(csvline{1}, RobName) % Suche Roboternamen (erste Spalte)
        % Zu modifizierenden Roboter gefunden. Diese Zeile verändert in Dateikopie schreiben
        found = true;
        % Zeile modifizieren
        csvline_mod = csvline;
        csvline_mod(1+8*N+3+(1:9)) = csvline_EE_FG_calc;
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
    fprintf('\tWert EE-FG auf %s gesetzt.\n', dec2bin(serroblib_csvline2bits_EE(csvline_EE_FG_calc),9));
  end
end
