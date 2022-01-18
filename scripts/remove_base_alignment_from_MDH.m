% Entferne die Basis-Ausrichtung aus der Datenbank. Zukünftig soll die
% Eigenschaft separat behandelt werden
% 
% Siehe auch: correct_phi_N_E.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

roblibpath=fileparts(which('serroblib_path_init.m'));
% Aktualisiere die Datenbank. Notwendig, damit nicht die gleichen
% Kinematiken nach der Änderung nochmal korrigiert werden.
serroblib_gen_bitarrays();
% Alle Roboter aus Datenbank laden
mdllistfile_Ndof = fullfile(roblibpath, 'serrob_list.mat');
l = load(mdllistfile_Ndof);
% Speichere die gefundenen Strukturen separat in einer Textdatei. 
outfile = fullfile(roblibpath, ['base_alignment_changed_list', ...
  datestr(now,'yyyymmdd_HHMMSS'), '.txt']);
fid_stats = fopen(outfile, 'w');
for j = 1:length(l.Names)
  RobName = l.Names{j};
  N = l.N(j);
  % Debug: Nur für bestimmte Kinematiken durchführen
  % if l.N(j) ~= 5, continue; end
  % if ~strcmp(RobName, 'S5RRRRR5'), continue; end
  RS = serroblib_create_robot_class(RobName);
  % Prüfe, ob MDH-Parameter zu Schema gehören
  RS.update_mdh(zeros(length(RS.pkin), 1));
  Tc = RS.fkine(zeros(RS.NJ,1));
  T1 = Tc(:,:,2);
  basedirt_x = false;
  basedirt_y = false;
  basedirt_z = false;
  if all(abs(abs(T1(1:3,3))-[1;0;0])<1e-10)
    basedirt_x = true;
  elseif all(abs(abs(T1(1:3,3))-[0;1;0])<1e-10)
    basedirt_y = true;
  elseif all(abs(abs(T1(1:3,3))-[0;0;1])<1e-10)
    basedirt_z = true;
  else
    error('Nicht vorgesehener Fall für Ausrichtung der Gelenkachse')
  end

  basedir_x = abs(RS.MDH.beta(1)) == pi/2 && abs(RS.MDH.alpha(1)) == pi/2;% && ...
%     (RS.MDH.sigma(1)==0 && abs(RS.MDH.offset(1)) == pi/2 || ...
%      RS.MDH.sigma(1)==1 && abs(RS.MDH.theta(1)) ==  pi/2);
  basedir_y = (RS.MDH.beta(1) == 0 && abs(RS.MDH.alpha(1)) == pi/2 || ...
              RS.MDH.beta(1) == pi && abs(RS.MDH.alpha(1)) == pi/2);% && ...
%     (RS.MDH.sigma(1)==0 && abs(RS.MDH.offset(1)) == pi/2 || ...
%      RS.MDH.sigma(1)==1 && abs(RS.MDH.theta(1)) ==  pi/2);
  basedir_z = RS.MDH.beta(1) == 0 && RS.MDH.alpha(1) == 0;
  if basedir_x ~= basedirt_x
    error('X-Basis-Ausrichtung stimmt nicht');
  end
  if basedir_y ~= basedirt_y
    error('Y-Basis-Ausrichtung stimmt nicht');
  end
  if basedir_z ~= basedirt_z
    error('Z-Basis-Ausrichtung stimmt nicht');
  end
  [csvline, csvbits] = serroblib_bits2csvline(l.BitArrays_Ndof(j,1));
  basedir_idx = find([basedirt_x;basedirt_y;basedirt_z]);
  fprintf('%d, %s: Base direction %s. joint %s, beta1=%s, alpha1=%s, theta1=%s, offset1=%s\n', ...
    j, RobName, char(119+basedir_idx), csvline{2}, csvline{3}, csvline{5}, ...
    csvline{7}, csvline{9});
   
  Mask_Origin = uint16(bin2dec('10000'));
  if bitand(Mask_Origin, l.BitArrays_Origin(j)) ~= 0
    fprintf('%d, %s: Manuell eingefügt. Überspringen.\n', j, RobName);
    continue
  end
  if basedir_z
    continue % Z-Ausrichtung ist Standard. Nichts unternehmen.
  end
  if strcmp(csvline{7}, 'theta1')
    fprintf('theta1 ist symbolisch\n');
  end

  %% Ändere die Werte in CSV-Tabelle
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
      csvline_mod{3} = '0'; % beta1 immer auf Null setzen
      csvline_mod{5} = '0'; % alpha1 auch
      if ~strcmp(csvline{7}, 'theta1')
        csvline_mod{7} = '0'; % theta1 auch, nur wenn konstant
      end
      csvline_mod{9} = '0'; % offset1 immer zu Null setzen. Hat keinen Einfluss
      % Aktualisiere auch die neue EE-Komponenten
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
  % Aktualisiere die mat-Datenbank zum Laden der aktualisierten Parameter.
  serroblib_gen_bitarrays(N);
  %% Erzeuge die symbolischen Funktionen neu
  serroblib_generate_mapleinput({RobName});
  serroblib_generate_code({RobName}, true, false, 1);
  fprintf(fid_stats, '%s\n', RobName);
  fprintf('Code für %d, %s neu generiert\n', j, RobName);
end
fclose(fid_stats);
% Aktualisiere die Datenbank
serroblib_gen_bitarrays();