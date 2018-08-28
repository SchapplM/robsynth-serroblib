% Generiere .mat-Dateien mit den Arrays zum schnelleren Durchsuchen der
% Datenbank

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function serroblib_gen_bitarrays(N_update)


if nargin == 0
  N_update = 1:6; % Aktualisiere alle Roboter
end
repopath=fileparts(which('serroblib_path_init.m'));
for N = N_update(:)'
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  BitArrays_Ndof = uint16(zeros(1,N));
  Names_Ndof = {};
  b = 1;
  % Durchsuche alle csv-Dateien im Ordner nach passenden Strukturen
  mdldir = fullfile(repopath, sprintf('mdl_%ddof', N));
  for d = dir(fullfile(mdldir, '*.csv'))'
    csvfilepath = fullfile(mdldir, d.name);
    fid = fopen(csvfilepath);
    tline = fgetl(fid);
    while ischar(tline)
      % Spaltenweise als Cell-Array
      csvline = strsplit(tline, ',');
      tline = fgetl(fid); % nächste Zeile
      if isempty(csvline) || strcmp(csvline{1}, '')
        continue
      end
      if length(csvline) < N*8+1
        % warning('Zeile %s sieht ungültig aus', tline_it);
        continue % nicht genug Spalten: Ungültiger Datensatz
      end
      % Prüfe, ob Zeile Roboterzeile ist
      robtype = d.name(1:end-4); % Name der csv-Datei ohne die Endung
      firstcol = csvline{1};
      if length(firstcol)<length(robtype) || ~strcmp(firstcol(1:length(robtype)), robtype)
        continue % keine Roboter-Zeile; wahrscheinlich Überschrift
      end
      BA = serroblib_csvline2bits(csvline);
      % dec2bin(BA)
      % fprintf('serroblib_gen_bitarrays: %s aus %s zur mat-Datei %s hinzugefügt\n', csvline{1}, d.name, sprintf('S%d_list.mat',N));
      % csvline2 = serroblib_bits2csvline(BA);

      Names_Ndof{b} = csvline{1}; %#ok<AGROW>
      BitArrays_Ndof(b,:) = BA;
      b = b+1;
      
    end
    fclose(fid);
  end
  
  % Alle Modelle neu in Ergebnisdatei (.mat) speichern
  % fprintf('serroblib_gen_bitarrays: Datei %s mit %d Einträgen gespeichert \n', mdllistfile_Ndof, size(BitArrays_Ndof,1));
  mkdirs(fileparts(mdllistfile_Ndof));
  save(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof');
end

