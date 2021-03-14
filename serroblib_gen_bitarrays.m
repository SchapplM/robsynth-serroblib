% Generiere .mat-Dateien mit den Arrays zum schnelleren Durchsuchen der
% Datenbank
% 
% Eingabe:
% N_update (optional)
%   Anzahl der Gelenk-Freiheitsgrade der Roboter, für die die .mat-Datei
%   aktualisiert werden soll.
% 
% Schreibt Dateien:
% mdl_xdof/Sx_list.mat (x=Anzahl FG aus Eingabe) mit Variablen:
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
%     4: Matlab-Code für dieses Modell liegt vor (1: Code für genau dieses
%        Modell liegt vor;0 : Kein Code vorliegend; 2: Code für
%        Haupt-Modell einer Variante liegt vor
%     5: Anzahl Drehgelenke. Damit kann schneller auf die prinzipielle
%        Struktur geschlossen werden (ist aber auch im Namen enthalten)
%     6: Anzahl technischer Gelenke. Ein mögliches Kugelgelenk zählt als 3,
%        ein Kardangelenk als 2, einfache Dreh- und Schubgelenke als 1
%     7: Reihenfolge technischer Gelenke. Dezimal kodiert.
%        1=R,2=P,3=C,4=U,5=S. 10er-Potenz entspricht Gelenknummer
%   BitArrays_Ndof_VF
%     Bit-Array mit Filter, um Varianten eines Roboters zu erkennen
%   BitArrays_Origin
%     Bit-Array mit Herkunft der Kinematik kodiert als Bits
% serrob_list.mat (wie mdl_xdof/Sx_list.mat, aber für alle Roboter)
%   Names
%     Namen der Roboter (anstatt Names_Ndof)
%   N
%     Anzahl der Gelenke des jeweiligen Roboters

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für Mechatronische Systeme, Universität Hannover

function serroblib_gen_bitarrays(N_update)


if nargin == 0
  N_update = 1:7; % Aktualisiere alle Roboter
end
repopath=fileparts(which('serroblib_path_init.m'));
mdllistfile_alldof = fullfile(repopath, 'serrob_list.mat');
for N = N_update(:)'
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  BitArrays_Ndof = uint16(zeros(0,N));
  BitArrays_Ndof_VF = uint16(zeros(0,N));
  BitArrays_phiNE = uint16(zeros(0,1));
  BitArrays_EEdof0 = uint16(zeros(0,1));
  BitArrays_Origin = uint16(zeros(0,1));
  AdditionalInfo = zeros(0,7);
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
      if length(csvline) ~= 1+N*8+3+6+3+2+5 % Siehe serroblib_gen_bitarrays.m
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
      [BAJ, BAR, BAE, BAJVF, BAO] = serroblib_csvline2bits(csvline);

      % Debug: Prüfe, die Hin- und Rückrechnung mit den Bit-Funktionen
%       [csvline_test] = serroblib_bits2csvline(BAJ);
%       I_nomatch = ~strcmp(csvline_test(2:end), csvline(2:length(csvline_test)));
%       if any(I_nomatch)
%         warning('%s: Die Hin- und Rückrechnung von CSV-Zeile und Bit-Arrays funktioniert nicht', Name);
%       end
      %% Weitere Eigenschaften bestimmen bzw. Auslesen
      % Lese die Nummer des letzten Positionsbestimmenden Gelenks aus
      lastposjoint = str2double(csvline{end-6});
      
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
        if isempty(variantof)
          warning('Modellvariante %s kann keinem Hauptmodell zugeordnet werden', Name);
          variantof = 0;
        end
      end
      
      % Prüfe, ob Matlab-Code vorliegt. Damit kann geprüft werden, ob eine
      % Variante einen eigenen Code-Ordner hat, oder den Code des
      % allgemeinen Haupt-Modells verwendet
      if isvariant
        code_dir_var = fullfile(repopath, sprintf('mdl_%ddof', N), name_haupt, ...
                                sprintf('hd_V%s', tokens_var{1}{4}));
        code_dir_gen = fullfile(repopath, sprintf('mdl_%ddof', N), name_haupt, ...
                                'hd');
        % Wenn nur Kinematik generiert wurde, liegen 16 m-Dateien vor.
        if length(dir(fullfile(code_dir_var, '*.m'))) > 20
          % Variante hat ihren eigenen Code-Ordner
          hascode = 1;
        elseif length(dir(fullfile(code_dir_gen, '*.m'))) > 20
          % Code für Variante nicht verfügbar. Haupt-Modell hat Code
          hascode = 2;
        else
          % Kein Code für Hauptmodell oder Variante
          hascode = 0;
        end
      else % Keine Variante
        code_dir = fullfile(repopath, sprintf('mdl_%ddof', N), Name, 'hd');
        % Wenn nur Kinematik generiert wurde, liegen 16 m-Dateien vor.
        if length(dir(fullfile(code_dir, '*.m'))) > 20
          % Code verfügbar (Genug Matlab-Funktionen im Ordner
          hascode = 1;
        else
          hascode = 0;
        end
      end
      
      % Anzahl Drehgelenke zählen
      numrotjoints = sum(Name=='R');
      
      % Anzahl technischer Gelenke bestimmen
      joints_string = csvline{end-5}; % die letzten 5 Spalten geben die Modellherkunft an
      if isempty(joints_string)
        % Ohne Angabe werden einzelne Dreh- und Schubgelenke angenommen
        % (stehen schon direkt im Namen drin)
        joints_string = Name(3:3+N-1);
      end
      numtechjoints = length(joints_string);
      % Technische Gelenke als Zahl kodieren
      factors = 10.^(0:numtechjoints-1);
      joints_number = 0;
      joints_number = joints_number + sum(1*double(joints_string=='R').*factors);
      joints_number = joints_number + sum(2*double(joints_string=='P').*factors);
      joints_number = joints_number + sum(3*double(joints_string=='C').*factors);
      joints_number = joints_number + sum(4*double(joints_string=='U').*factors);
      joints_number = joints_number + sum(5*double(joints_string=='S').*factors);
      
      %% Ausgabe belegen
      Names_Ndof{b} = csvline{1}; %#ok<AGROW>
      BitArrays_Ndof(b,:) = BAJ;
      BitArrays_Ndof_VF(b,:) = BAJVF;
      BitArrays_phiNE(b,:) = BAR;
      BitArrays_EEdof0(b,:) = BAE;
      BitArrays_Origin(b,:) = BAO;
      AdditionalInfo(b,:) = [lastposjoint, double(isvariant), variantof, hascode, numrotjoints, numtechjoints, joints_number];
      b = b+1;
    end
    fclose(fid);
  end
  % Alle Modelle neu in Ergebnisdatei (.mat) speichern
  Ndof_struct_new = struct( ...
    'Names_Ndof', {Names_Ndof}, 'BitArrays_Ndof',  BitArrays_Ndof, ...
    'BitArrays_Ndof_VF', BitArrays_Ndof_VF, 'BitArrays_phiNE',  BitArrays_phiNE, ...
    'BitArrays_EEdof0', BitArrays_EEdof0, 'BitArrays_Origin',  BitArrays_Origin, ...
    'AdditionalInfo', AdditionalInfo);
  % Prüfe, ob die Datei schon existiert und ob sie sich ändern würde. Ver-
  % meide das Schreiben identischer Daten in die Datei. Durch parallele
  % Schreibprozesse kann die Datei beschädigt werden (Rechencluster)
  write_new = false;
  if exist(mdllistfile_Ndof, 'file')
    Ndof_struct_old = load(mdllistfile_Ndof);
    if ~isequaln(Ndof_struct_old, Ndof_struct_new)
      write_new = true; % Dateiinhalt wird sich ändern. Schreibe neu.
    end
  else
    mkdirs(fileparts(mdllistfile_Ndof)); % Ordner könnte nicht existieren
    write_new = true; % Schreibe mat-Datei neu, sie existiert noch nicht
  end
  if write_new
    save(mdllistfile_Ndof, '-struct', 'Ndof_struct_new');
%     fprintf('serroblib_gen_bitarrays: Datei %s mit %d Einträgen gespeichert \n', ...
%       mdllistfile_Ndof, size(BitArrays_Ndof,1));
  end
  % Gesamt-Datei aktualisieren
  if length(N_update) == 7 % nur wenn auch alle Roboter aktualisiert werden
    % Anzahl der Spalten erhöhen (Ndof-Variable wird nicht mehr benötigt
    Ndof_struct_new.BitArrays_Ndof = [Ndof_struct_new.BitArrays_Ndof, ...
      uint8(zeros(size(Ndof_struct_new.BitArrays_Ndof,1), 7-N))];
    Ndof_struct_new.BitArrays_Ndof_VF = [Ndof_struct_new.BitArrays_Ndof_VF, ...
      uint8(zeros(size(Ndof_struct_new.BitArrays_Ndof_VF,1), 7-N))];
    if N == 1
      alldof_struct_new = Ndof_struct_new;
      alldof_struct_new.N = ones(length(Ndof_struct_new.Names_Ndof),1);
      alldof_struct_new.Names = alldof_struct_new.Names_Ndof;
      alldof_struct_new = rmfield(alldof_struct_new, 'Names_Ndof');
    else
      alldof_struct_new.Names = [alldof_struct_new.Names, ...
        Ndof_struct_new.Names_Ndof];
      for e = {'BitArrays_Ndof', 'BitArrays_Ndof_VF', 'BitArrays_phiNE', ...
          'BitArrays_EEdof0', 'BitArrays_Origin', 'AdditionalInfo'}
        alldof_struct_new.(e{1}) = [alldof_struct_new.(e{1}); Ndof_struct_new.(e{1})];
      end
      alldof_struct_new.N = [alldof_struct_new.N; repmat(N,length(Ndof_struct_new.Names_Ndof),1)];
    end
  end
end
% Gesamt-Datei speichern (falls vollständig zusammengestellt wurde)
if length(N_update) == 7
  write_new = false;
  if exist(mdllistfile_alldof, 'file')
    alldof_struct_old = load(mdllistfile_alldof);
    if ~isequaln(alldof_struct_old, alldof_struct_new)
      write_new = true; % Dateiinhalt wird sich ändern. Schreibe neu.
    end
  else
    write_new = true; % Schreibe mat-Datei neu, sie existiert noch nicht
  end
  if write_new
    save(mdllistfile_alldof, '-struct', 'alldof_struct_new');
  end
end