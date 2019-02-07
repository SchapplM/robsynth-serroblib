% Entferne ein Robotermodell aus der Datenbank (csv-Tabelle)
% Die Zeile des Robotermodells (z.B. RRPR4) wird aus der csv-Tabelle der
% Gelenkstruktur (z.B. RRPR) entfernt. Zusätzlich werden die
% Matlab-Funktionen entfernt.
% 
% Eingabe:
% mdlname
%   Name des Roboters in der Datenbank
% 
% Ausgabe:
% success
%   true, falls Roboter erfolgreich entfernt wurde

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

function success = serroblib_remove_robot(Name)

% Roboterdaten aus Namen extrahieren
N = str2double(Name(2));
typestring = Name(3:3+N-1);

% Dateinamen der csv-Tabellen bestimmen
filename = sprintf('S%d%s.csv', N, typestring);
repopath=fileparts(which('serroblib_path_init.m'));
filepath_csv = fullfile(repopath, sprintf('mdl_%ddof', N), filename);
filepath_csv_copy = [filepath_csv, '.copy']; % Kopie der Tabelle zur Bearbeitung
% CSV-Tabelle zeilenweise durchgehen
fid = fopen(filepath_csv);
fidc = fopen(filepath_csv_copy, 'w');
tline = fgetl(fid);
found = false;
while ischar(tline)
  % Spaltenweise als Cell-Array
  csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
  if strcmp(csvline{1}, Name)
    % Zu löschenden Roboter gefunden. Diese Zeile nicht in Dateikopie schreiben
    found = true;
  else
    % Zeile in die Dateikopie schreiben
    fwrite(fidc, [tline, newline()]);
  end
  tline = fgetl(fid); % nächste Zeile
end
fclose(fid);
fclose(fidc);
if ~found
  success = false;
  warning('Zu löschendes Modell %s nicht in %s gefunden', ...
    Name, filepath_csv);
  return
end
% Modifizierte Tabelle zurückkopieren
copyfile(filepath_csv_copy, filepath_csv);
% Kopie-Tabelle löschen
delete(filepath_csv_copy);
% Ordner mit Code und Parameter-Modellen löschen
% Aktuell muss noch manuell geprüft werden, ob dabei wichtige Daten verloren gehen
robdir = fullfile(repopath, sprintf('mdl_%ddof', N), Name);
rmdir(robdir, 's');

success = true;