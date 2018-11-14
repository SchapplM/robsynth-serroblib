% Erstelle csv-Tabellen für eine mit Namen gegebene Roboterstruktur
% (als Struktur wird die Gelenkanordnung (PRRR) mit symbolisch definierten
% MDH-Parametern verstanden)
% In die Tabelle sollen dann alle bekannten Roboterstrukturen mit
% Zahlenwerten für die MDH-Parameter eingetragen werden
% 
% Eingabe:
% Name
%   Name des Roboter, für die die csv-Tabelle erstellt werden soll
%   Beispiel: 'S3RRR1'
% 
% Ausgabe:
% exists
%   true, falls die csv-Datei schon existiert (dann keine Aktion)
% filepath_csv
%   Pfad der zu erstellenden csv-Datei
% 
% Schreibt Dateien:
% /mdl_xdof/SRR...PRy/models.csv (für den gegebenen Roboter)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function [exists, filepath_csv] = serroblib_create_robot_csv(Name)

exists = false;

N = str2double(Name(2)); % Schema: SxRRRR4; x="Anzahl FG"
% Typ = Name(1:N+2); % Nehme Teil SxRRRR ohne die laufende Nummer
repopath=fileparts(which('serroblib_path_init.m'));
filepath_csv = fullfile(repopath, sprintf('mdl_%ddof', N), Name, sprintf('models.csv'));

% Wenn Datei schon existiert, erstelle nicht neu
if exist(filepath_csv, 'file')
  exists = true;
  return
end

% Parameter des Roboters laden
mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof');
BA = l.BitArrays_Ndof(strcmp(l.Names_Ndof, Name),:);
csvline = serroblib_bits2csvline(BA);
%% Erstelle Überschriftenzeilen für die Struktur
% Siehe auch: serroblib_add_robot.m
% Die ersten drei Spalten sind zur Kennzeichnung des Roboters
% Erste Zeile Gesamtüberschrift
csvline_head1 = {'Name/lfdNr', 'Name2', 'Langname'};
% Zweite Zeile: Bezeichnung der MDH-Parameter
csvline_head2 = {'', '', ''};
% Dritte Zeile: Kennzeichnung, welche MDH-Parameter Null oder fest auf z.B.
% pi/2 liegen
csvline_head3 = {Name, '', ''};
c = 4; % Spalten-Zähler
% Kopfzeile für alle Gelenke erstellen
for kk = 1:N
  csvline_head1{c} = sprintf('Gelenk %d', kk);
  for jj = 1:12
    csvline_head1{c+jj} = '';
  end
  for jj = 1:12
    if jj < 9
      csvline_head3{c+jj-1} =  csvline{1+(kk-1)*8+jj};
    else
      csvline_head3{c+jj-1} = '';
    end
  end
  csvline_head2{c} = 'Typ'; c = c+1;
  csvline_head2{c} = 'beta'; c = c+1;
  csvline_head2{c} = 'b'; c = c+1;
  csvline_head2{c} = 'alpha'; c = c+1;
  csvline_head2{c} = 'a'; c = c+1;
  csvline_head2{c} = 'theta'; c = c+1;
  csvline_head2{c} = 'd'; c = c+1;
  csvline_head2{c} = 'offset'; c = c+1;
  csvline_head2{c} = 'qmin'; c = c+1;
  csvline_head2{c} = 'qmax'; c = c+1;
  csvline_head2{c} = 'vmax'; c = c+1;
  csvline_head2{c} = 'VZ'; c = c+1;
  csvline_head2{c} = 'just'; c = c+1;
end

% Kopfzeile für zusätzliche Infos
c=c+1;
csvline_head1{c} = 'Allgemein';
csvline_head2{c} = 'Nenn-Traglast';
csvline_head3{c} = '';
c=c+1;
csvline_head1{c} = '';
csvline_head2{c} = 'Masse';
csvline_head3{c} = '';

% String aus Cell-Arrays erzeugen
line_head1 = csvline_head1{1};
line_head2 = csvline_head2{1};
line_head3 = csvline_head3{1};
for i = 2:length(csvline_head1)
  line_head1 = sprintf('%s;%s', line_head1, csvline_head1{i});
  line_head2 = sprintf('%s;%s', line_head2, csvline_head2{i});
  line_head3 = sprintf('%s;%s', line_head3, csvline_head3{i});
end
% Kopfzeilen in csv-Tabelle schreiben
mkdirs(fileparts(filepath_csv));
fid = fopen(filepath_csv, 'a');
fwrite(fid, [line_head1, newline]);
fwrite(fid, [line_head2, newline]);
fwrite(fid, [line_head3, newline]);
fclose(fid);