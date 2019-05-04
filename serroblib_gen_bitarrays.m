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
%   BitArrays_phiNE
%     Bit-Arrays mit Transformation KS-N zu KS-E (in Euler-Winkeln)
%     Diese Information dient nur zur Konsistenz der Kinematik
%   AdditionalInfo
%     Vektor mit zusätzlichen Eigenschaften des Roboters (als Zahl).
%     Einträge:
%     1: Letztes positionsbeeinflussendes Gelenk
%     2: Ist Variante (1/0)
%     3: Nummer des Modells, von dem dieses Modell eine Variante ist
%        Die Nummer bezieht sich auf die Liste aller Modelle und Varianten
%        in der Tabelle
%   BitArrays_Ndof_VF
%     Bit-Array mit Filter, um Varianten eines Roboters zu erkennen

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für Mechatronische Systeme, Universität Hannover

function serroblib_gen_bitarrays(N_update)


if nargin == 0
  N_update = 1:7; % Aktualisiere alle Roboter
end
repopath=fileparts(which('serroblib_path_init.m'));
for N = N_update(:)'
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  BitArrays_Ndof = uint16(zeros(1,N));
  BitArrays_Ndof_VF = uint16(zeros(1,N));
  BitArrays_phiNE = uint16(zeros(1,1));
  BitArrays_EEdof0 = uint16(zeros(1,1));
  AdditionalInfo = zeros(1,3);
  Names_Ndof = {};
  b = 1; % Zähler für gefundene Roboterkonfigurationen aus csv-Tabelle für N FG
  %% Durchsuche alle csv-Dateien im Ordner nach passenden Strukturen
  mdldir = fullfile(repopath, sprintf('mdl_%ddof', N));
  for d = dir(fullfile(mdldir, '*.csv'))'
    robtype = d.name(1:end-4); % Name der csv-Datei ohne die Endung
    csvfilepath = fullfile(mdldir, d.name);
    fid = fopen(csvfilepath);
    tline = fgetl(fid);
    while ischar(tline)
      %% CSV-Zeile verarbeitejn und filtern
      % Spaltenweise als Cell-Array
      csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
      tline = fgetl(fid); % nächste Zeile
      if isempty(csvline) || strcmp(csvline{1}, '')
        continue
      end
      if length(csvline) ~= 1+N*8+3+6+3+1
        warning('Zeile %s (Datei %s) sieht ungültig aus', tline, d.name);
        continue % nicht genug Spalten: Ungültiger Datensatz
      end
      % Prüfe, ob Zeile Roboterzeile ist
      firstcol = csvline{1};
      if length(firstcol)<length(robtype) || ~strcmp(firstcol(1:length(robtype)), robtype)
        continue % keine Roboter-Zeile; wahrscheinlich Überschrift
      end
      % Ab hier ist davon auszugehen, dass die Zeile eine Roboter-Zeile ist
      Name = csvline{1};
      
      % Roboterzeile in Binär-kodierte Informationen umrechnen
      [BAJ, BAR, BAE, BAJVF] = serroblib_csvline2bits(csvline);

      %% Weitere Eigenschaften bestimmen bzw. Auslesen
      % Lese die Nummer des letzten Positionsbestimmenden Gelenks aus
      lastposjoint = str2double(csvline{end});
      
      % Bestimme, ob der Roboter eine Variante oder ein Hauptmodell ist
      exp_mdl = 'S(\d)([RP]+)(\d+)'; % Format "S3RRR1" für Hauptmodelle
      exp_var = 'S(\d)([RP]+)(\d+)V(\d+)'; % Format "S3RRR1V1" für Varianten
      [tokens_mdl, ~] = regexp(Name,exp_mdl,'tokens','match');
      [tokens_var, ~] = regexp(Name,exp_var,'tokens','match');
      if ~isempty(tokens_mdl) && isempty(tokens_var) % Name entspricht Haupt-Modell
        isvariant = false;
      elseif ~isempty(tokens_var) % Name entspricht Variante
        isvariant = true;
      else
        error('Namensschema %s passt weder zu Hauptmodellen noch zu Variantenmodellen', Name);
      end
      
      % Bestimme die Nummer des Roboters in der Liste, von dem die Variante
      % abstammt
      if ~isvariant
        variantof = b;
      else
        name_haupt = ['S',tokens_var{1}{1},tokens_var{1}{2},tokens_var{1}{3}];
        variantof = find(strcmp(Names_Ndof, name_haupt));
      end
      %% Ausgabe belegen
      Names_Ndof{b} = csvline{1}; %#ok<AGROW>
      BitArrays_Ndof(b,:) = BAJ;
      BitArrays_Ndof_VF(b,:) = BAJVF;
      BitArrays_phiNE(b,:) = BAR;
      BitArrays_EEdof0(b,:) = BAE;
      AdditionalInfo(b,:) = [lastposjoint, double(isvariant), variantof];
      b = b+1;
    end
    fclose(fid);
  end
  
  % Alle Modelle neu in Ergebnisdatei (.mat) speichern
  % fprintf('serroblib_gen_bitarrays: Datei %s mit %d Einträgen gespeichert \n', mdllistfile_Ndof, size(BitArrays_Ndof,1));
  mkdirs(fileparts(mdllistfile_Ndof));
  save(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_Ndof_VF', ...
    'BitArrays_phiNE', 'BitArrays_EEdof0', 'AdditionalInfo');
end
