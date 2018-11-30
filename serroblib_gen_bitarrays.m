% Generiere .mat-Dateien mit den Arrays zum schnelleren Durchsuchen der
% Datenbank
% 
% Eingabe:
% N_update (optional)
%   Anzahl der Gelenk-Freiheitsgrade der Roboter, für die die .mat-Datei
%   aktualisiert werden soll.
% 
% Schreibt Dateien:
% mdl_xdof.mat (x=Anzahl FG aus Eingabe) mit Variablen:
%   Names_Ndof
%     Cell-Array mit Namen aller Roboter dieses FG
%   BitArrays_Ndof
%     Bit-Arrays mit gespeicherten Eigenschaften der einzelnen Gelenke
%   BitArrays_EEdof0
%     Bit-Arrays mit Eigenschaften der Endeffektor-Bewegung (im Basis-KS)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function serroblib_gen_bitarrays(N_update)


if nargin == 0
  N_update = 1:7; % Aktualisiere alle Roboter
end
repopath=fileparts(which('serroblib_path_init.m'));
for N = N_update(:)'
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  BitArrays_Ndof = uint16(zeros(1,N));
  BitArrays_EEdof0 = uint16(zeros(1,1));
  Names_Ndof = {};
  b = 1; % Zähler für gefundene Roboterkonfigurationen aus csv-Tabelle für N FG
  % Durchsuche alle csv-Dateien im Ordner nach passenden Strukturen
  mdldir = fullfile(repopath, sprintf('mdl_%ddof', N));
  for d = dir(fullfile(mdldir, '*.csv'))'
    % fprintf('%s\n', d.name);
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
        % warning('Zeile %s sieht ungültig aus', tline_it);
        continue % nicht genug Spalten: Ungültiger Datensatz
      end
      % Prüfe, ob Zeile Roboterzeile ist
      robtype = d.name(1:end-4); % Name der csv-Datei ohne die Endung
      firstcol = csvline{1};
      if length(firstcol)<length(robtype) || ~strcmp(firstcol(1:length(robtype)), robtype)
        continue % keine Roboter-Zeile; wahrscheinlich Überschrift
      end
%       % Entferne Spalten, die keinen Gelenk-Bezug haben
%       csvline_jointkin = csvline(1:end-6);
      [BAJ, BAE] = serroblib_csvline2bits(csvline);

      % Ausgabe belegen
      Names_Ndof{b} = csvline{1}; %#ok<AGROW>
      BitArrays_Ndof(b,:) = BAJ;
      BitArrays_EEdof0(b,:) = BAE;
      b = b+1;
    end
    fclose(fid);
  end
  
  % Alle Modelle neu in Ergebnisdatei (.mat) speichern
  % fprintf('serroblib_gen_bitarrays: Datei %s mit %d Einträgen gespeichert \n', mdllistfile_Ndof, size(BitArrays_Ndof,1));
  mkdirs(fileparts(mdllistfile_Ndof));
  save(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_EEdof0');
end

