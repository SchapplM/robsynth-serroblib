% Finde ein Robotermodell in der Datenbank
% 
% Eingabe:
% csvline {1xM} cell array
%   Cell array mit den Spalten der csv-Tabelle zu dem gesuchten Roboter
%   Modellname
%   Spalten für alle MDH-Parameter. Siehe serroblib_add_robot.m
% filter_relkin
%   Aktivieren des Filters für die reine Betrachtung der relativen
%   Kinematik zwischen den Gelenken. Beim Finden des Roboters ist es egal,
%   ob die erste Achse in z-Richtung verläuft oder gedreht ist.
% filter_var
%   prüfe, ob neue Struktur eine Variante eines bestehenden Hauptmodells
%   ist
% 
% Rückgabe:
% found
%   true, wenn gegebener Roboter bereits in der .mat-Datenbank vorliegt
% idx [4x1] (index_all, index_type, index_gen, index_var)
%   1: Nummer des Roboters in der Liste aller Roboter mit dieser
%      Gelenkanzahl
%   2: Nummer des Roboters in der Liste aller Roboter dieser Gelenkreihenfolge
%      (falls gefunden). Nummer des letzten freien Roboters, falls nicht gefunden.
%      Die Liste beinhaltet Hauptmodelle und Varianten
%   3: Nummer des gefundenen Roboters in allen Hauptmodellen. Falls nicht
%      gefunden: Erste freie Nummer.
%   4: Nummer des gefundenen Roboters in den Varianten des Hauptmodells. Falls
%      nicht gefunden: Anzahl der Varianten dieses Hauptmodells
% num [2x1] (num_all, num_gen)
%   1: Anzahl der Roboter dieser Gelenkreihenfolge insgesamt (bezogen auf
%      Hauptmodelle und Varianten)
%   2: Anzahl der Roboter in der Liste der Hauptmodelle dieser
%      Gelenkreihenfolge (ohne Varianten)
% Name
%   Name des durch Parameter gegebenen Roboters (falls gefunden)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [found, idx, num, Name] = serroblib_find_robot(csvline, filter_relkin, filter_var)

%% Initialisierung
found = false;
index_all = 0;
index_type = 0;
index_gen = 0;
index_var = 0;
Name = '';
% Reguläre Ausdrücke zur Erkennung von Modellen und Varianten. Benutze
% Anfangs- und End-Begrenzung ("^","$"), damit die Variante nicht als
% Hauptmodell erkannt wird
exp_mdl = '^S(\d)([RP]+)(\d+)$'; % Format "S3RRR1" für Hauptmodelle
exp_var = '^S(\d)([RP]+)(\d+)V(\d+)$'; % Format "S3RRR1V1" für Varianten

if nargin < 2
  filter_relkin = true;
end
if nargin < 3
  filter_var = false;
end

% csv-Zeile in Bit-Array umwandeln
BA = serroblib_csvline2bits(csvline);

%% Bit-Array in Liste suchen
N = length(BA); % Anzahl Gelenke

% Daten laden
repopath=fileparts(which('serroblib_path_init.m'));
mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_Ndof_VF');

% Filter-Variable für erstes Gelenk: Die ersten Parameter beta, b, alpha,
% a, d, theta sollen nicht ausgewertet werden
BAJ_Iso_Filter = uint16(zeros(1,N));
for i = 1:N
  BAJ_Iso_Filter(i) = intmax;
end
if filter_relkin
  % Setze alle Bits auf 1, bis auf die Bits der Orientierung der ersten Achse
  % (Zu den Bits, siehe serroblib_bits2csvline.m)
  for i = 2:13 % Bits von beta bis theta
    BAJ_Iso_Filter(1) = bitset(BAJ_Iso_Filter(1), i, 0);
  end % dec2bin(BAJ_Iso_Filter)
  % Bit für Offset Null setzen (für alle Gelenke)
  for aa = 1:N
    % Offset Filtern: Der Offset auf die Gelenkkoordinate q ist egal für die
    % Kinematik (im Sinne von Isomorphismen)
    for i = 14:16
      BAJ_Iso_Filter(aa) = bitset(BAJ_Iso_Filter(aa), i, 0);
    end
  end
else
  % Filter ist nicht aktiv, da alle Bits des Filters auf "1" gelassen
  % werden.
end

% Filter für Typ des Roboters und DH-Parameter
num_all = 0;
num_gen = 0;
num_var = 0;
firstfree_gen = 0;
idx_gen_last = 0;
for i = 1:size(l.BitArrays_Ndof, 1) % Alle Roboterstrukturen aus Datenbank durchgehen
  % Zähle, wie viele Roboter dieser Gelenkreihenfolge existieren
  % (das erste Bit kennzeichnet den Gelenktyp)
  if all(bitand(BA, 1) == bitand(l.BitArrays_Ndof(i,:), 1))
    num_all = num_all + 1;
    % Zähle, wie viele Haupt-Modelle das waren
    [tokens_mdl, ~] = regexp(l.Names_Ndof{i},exp_mdl,'tokens','match');
    if ~isempty(tokens_mdl)
      num_gen = num_gen + 1;
      index_var = 0; % Neue Hauptstruktur führt zum Zurücksetzen des Zählers für Varianten
      % Prüfe, ob die fortlaufende Nummerierung stimmt
      if firstfree_gen == 0 && str2double(tokens_mdl{1}{3}) ~= idx_gen_last + 1
        firstfree_gen = idx_gen_last + 1;
      end
      idx_gen_last = str2double(tokens_mdl{1}{3}); % Nummer des vorherigen Modells für nächsten Durchlauf
    end
    % Zähle, wie viele Varianten dieses Haupt-Modells es waren
    % Annahme: Hauptmodelle und Varianten stehen alle untereinander
    [tokens_var, ~] = regexp(l.Names_Ndof{i},exp_var,'tokens','match');
    if ~isempty(tokens_var)
      num_var = num_var + 1;
      index_var = str2double(tokens_var{1}{4});
    end
  end
  
  % Prüfe ob alle Bits übereinstimmen (alle MDH-Parameter)
  % Nutze dazu den Filter für das erste Gelenk
  if(all( bitand(BA,BAJ_Iso_Filter) == bitand(l.BitArrays_Ndof(i,:),BAJ_Iso_Filter) ))
    % Eintrag ist die gesuchte Nummer (bezogen auf die Roboter derselben
    % Gelenkreihenfolge)
    index_type = num_all; % Index in den Robotern mit derselben Gelenkreihenfolge
    index_gen = num_gen; % Index in den Hauptmodellen mit derselben Gelenkreihenfolge
    index_all = i; % Index in allen Robotern mit derselben Anzahl Gelenk-FG
    found = true;
    break
  end
end
if ~found
  if firstfree_gen == 0
    % Es wurde keine Lücke zum Eintragen gefunden. Die nächste freie Nummer
    % wird ans Ende angehängt
    firstfree_gen = num_gen + 1;
  end
  % Rückgabe des letzten freien Index dieser Gelenkreihenfolge
  index_gen = firstfree_gen;
end
%% Variante suchen

if ~found && filter_var
  zlr_type = 0; % Zähler für die Position der Variante in allen Robotern (dieser FG-Zahl)
  zlr_gen = 0; % Zähler für Nummer des Allgemeinen Modells der Variante
  index_var = 0;
  for i = 1:size(l.BitArrays_Ndof, 1) % Alle Roboterstrukturen aus Datenbank durchgehen
    if ~all(bitand(BA, 1) == bitand(l.BitArrays_Ndof(i,:), 1))
      continue % Gelenkreihenfolge stimmt nicht. Weitergehen
    end
    zlr_type = zlr_type + 1; % Gleiche Gelenkfolge: Hochzählen (entspricht Nummer in Datei
    % Erster Schritt: Suche, welches Hauptmodell es ist
    BAJ_Filter_Ges = bitand(BAJ_Iso_Filter, l.BitArrays_Ndof_VF(i,:)); % Filter für Hauptmodell
    [tokens_mdl, ~] = regexp(l.Names_Ndof{i},exp_mdl,'tokens','match');
    if ~isempty(tokens_mdl) && ~found
      zlr_gen = zlr_gen + 1; % Ist Hauptmodell
    end
    if ~found && (all( bitand(BA,BAJ_Filter_Ges) == bitand(l.BitArrays_Ndof(i,:),BAJ_Filter_Ges) ))
      % Eintrag stimmt überein
      index_type = zlr_type; % Index in den Robotern mit derselben Gelenkreihenfolge
      index_all = i; % Index in allen Robotern mit derselben Anzahl Gelenk-FG
      found = true;
      [tokens_mdl_found, ~] = regexp(l.Names_Ndof{i},exp_mdl,'tokens','match');
      index_gen = str2double(tokens_mdl_found{1}{3});
      continue % sonst wird die Prüfung der folgenden Zeilen schon auf diese Zeile angewendet
    end
    % Zweiter Schritt: Finde heraus, wie viele Varianten dieses Hauptmodell
    % schon hat
    if found
      [tokens_var, ~] = regexp(l.Names_Ndof{i},exp_var,'tokens','match');
      if isempty(tokens_var) && ...
          (~strcmp(tokens_mdl{1}{2}, tokens_mdl_found{1}{2}) || ~strcmp(tokens_mdl{1}{3}, tokens_mdl_found{1}{3}))
        % Der aktuelle Roboter ist ein Hauptmodell, dass sich von dem         
        % gefundenen unterscheidet. Die Suche nach Varianten des Modells
        % kann abgebrochen werden
        break;
      end
      if ~isempty(tokens_var) && ...
          (~strcmp(tokens_mdl_found{1}{2}, tokens_var{1}{2}) || ~strcmp(tokens_mdl_found{1}{3}, tokens_var{1}{3}))
        % Aktueller Roboter ist Variante, aber keine Variante mehr des
        % gefundenen Hauptmodells. Höre auf zu suchen 
        break;
      end
      % Bestimme die Anzahl der Varianten, die es für den Roboter gibt
      index_var = str2double(tokens_var{1}{4});
    end
  end
end

%% Rückgabewerte belegen
if found
  Name = l.Names_Ndof{index_all};
end
idx = [index_all; index_type; index_gen; index_var];
num = [num_all; num_gen];