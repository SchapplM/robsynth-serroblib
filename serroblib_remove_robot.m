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
%   true, falls Roboter aus Datenbank entfernt wurde

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
else
  success = true;
  % Modifizierte Tabelle zurückkopieren
  copyfile(filepath_csv_copy, filepath_csv);
  % Kopie-Tabelle löschen
  delete(filepath_csv_copy);
end

%% Ordner mit Code und Parameter-Modellen löschen
% Aktuell muss noch manuell geprüft werden (über git), ob dabei wichtige
% Daten verloren gehen. Die Löschung wird vorgenommen, auch wenn der
% Roboter nicht in der Datenbank war.

% Prüfe, ob Modellname eine Variante ist
exp_var = '^S(\d)([RP]+)(\d+)V(\d+)$'; % Format "S3RRR1V1" für Varianten
[tokens_var, ~] = regexp(Name,exp_var,'tokens','match');

if isempty(tokens_var)
  % Keine Variante: Lösche komplettes Verzeichnis
  robdir = fullfile(repopath, sprintf('mdl_%ddof', N), Name);
  if exist(robdir, 'file')
    rmdir(robdir, 's');
    fprintf('Ordner %s entfernt\n', robdir);
  end
else
  % Variante: Lösche nur Code-Ordner der Variante
  Name_GenMdl = sprintf('S%s%s%s', tokens_var{1}{1}, tokens_var{1}{2}, tokens_var{1}{3});
  fcn_dir_var = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, ...
    sprintf('hd_V%s', tokens_var{1}{4}));
  if exist(fcn_dir_var, 'file')
    rmpath(fcn_dir_var);
    fprintf('Ordner %s entfernt\n', fcn_dir_var);
  end
  % Lösche außerdem Konvertierungsfunktionen, falls die Variante keinen
  % Code-Ordner hat
  fcnlist = {'pkin_gen2var', 'pkin_var2gen'};
  for f = fcnlist
    fcndat = fullfile(fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, ...
      'var', [Name, '_', f{1}, '.m']));
    if exist(fcndat, 'file')
      delete(fcndat);
    end
  end
end
