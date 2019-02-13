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
% 
% Rückgabe:
% found
%   true, wenn gegebener Roboter bereits in der .mat-Datenbank vorliegt
% index
%   Nummer des Roboters in der Liste aller Roboter dieser Gelenkreihenfolge
%   (falls gefunden) 
% num
%   Anzahl der Roboter dieser Gelenkreihenfolge insgesamt
% Name
%   Name des durch Parameter gegebenen Roboters (falls gefunden)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function [found, index, num, Name] = serroblib_find_robot(csvline, filter_relkin)

%% Initialisierung
found = false;
index = 0;
Name = '';

if nargin < 2
  filter_relkin = true;
end

% csv-Zeile in Bit-Array umwandeln
BA = serroblib_csvline2bits(csvline);

%% Bit-Array in Liste suchen
N = length(BA); % Anzahl Gelenke

% Daten laden
repopath=fileparts(which('serroblib_path_init.m'));
mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof');

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
num = 0;
for i = 1:size(l.BitArrays_Ndof, 1) % Alle Roboterstrukturen aus Datenbank durchgehen
  % Zähle, wie viele Roboter dieser Gelenkreihenfolge existieren
  % (das erste Bit kennzeichnet den Gelenktyp)
  num = num + all(bitand(BA, 1) == bitand(l.BitArrays_Ndof(i,:), 1));
  
  % Prüfe ob alle Bits übereinstimmen (alle MDH-Parameter)
  % Nutze dazu den Filter für das erste Gelenk
  if(all( bitand(BA,BAJ_Iso_Filter) == bitand(l.BitArrays_Ndof(i,:),BAJ_Iso_Filter) ))
    % Eintrag ist die gesuchte Nummer (bezogen auf die Roboter derselben
    % Gelenkreihenfolge)
    index_type = num; %Index in den Robotern mit derselben Gelenkreihenfolge
    index_N = i; % Index in allen Robotern mit derselben Anzahl Gelenk-FG
    found = true;
    break
  end
end
% Rückgabewerte belegen
if found
  index = index_type;
  Name = l.Names_Ndof{index_N};
end