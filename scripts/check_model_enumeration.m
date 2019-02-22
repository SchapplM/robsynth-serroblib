% Prüfe, ob alle Modelle in der Datenbank fortlaufend nummeriert sind.
% Durch Löschung kann es passieren, dass Lücken in der Nummerierung
% auftreten. Z.B. wenn in SPPRR.csv nach SPPRR1 direkt SPPRR3 kommt.
% 
% Ergebnis:
% Fehlermeldung, wenn Nummerierung nicht stimmt

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Initialisierung
roblibpath=fileparts(which('serroblib_path_init.m'));

%% Durchsuche alle csv-Tabellen und prüfe die Reihenfolge
zlr_ges = 0;
for N = 1:7
  zlr_ges_N = 0;
  mdldir = fullfile(roblibpath, sprintf('mdl_%ddof', N));
  for d = dir(fullfile(mdldir, '*.csv'))'
    lfdnr_vorher = 0;
    robtype = d.name(1:(N+2)); % Name der csv-Datei ohne die Endung
    csvfilepath = fullfile(mdldir, d.name);
    fid = fopen(csvfilepath);
    tline = fgetl(fid);
    while ischar(tline)
      % Spaltenweise als Cell-Array
      csvline = strsplit(tline, ';');
      tline = fgetl(fid); % nächste Zeile
      if isempty(csvline) || strcmp(csvline{1}, '')
        continue
      end
      if length(csvline) < N*8+1+6
        continue % nicht genug Spalten: Ungültiger Datensatz
      end
      % Prüfe, ob Zeile Roboterzeile ist
      firstcol = csvline{1};
      if length(firstcol)<length(robtype) || ~strcmp(firstcol(1:length(robtype)), robtype)
        continue % keine Roboter-Zeile; wahrscheinlich Überschrift
      end
      lfdnr=str2double(csvline{1}(length(robtype)+1:end));
      if lfdnr_vorher ~= lfdnr-1
        error('Laufende Nummern passen nicht: %s', firstcol);
      end
      % Belegung für nächste Iteration
      lfdnr_vorher = lfdnr;
      zlr_ges_N = zlr_ges_N + 1;
      zlr_ges = zlr_ges + 1;
    end
    fclose(fid);
  end
  fprintf('%d Roboterstrukturen mit %d FG geprüft. Reihenfolge stimmt.\n', zlr_ges_N, N);
end
fprintf('%d Roboterstrukturen insgesamt geprüft. Reihenfolge stimmt.\n', zlr_ges);