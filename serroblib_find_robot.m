% Finde ein Robotermodell in der Datenbank
% 
% Eingabe:
% csvline {1xM} cell array
%   Cell array mit den Spalten der csv-Tabelle zu dem gesuchten Roboter
%   Modellname
%   Spalten für alle MDH-Parameter. Siehe serroblib_add_robot.m
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

function [found, index, num, Name] = serroblib_find_robot(csvline)

%% Initialisierung
found = false;
index = 0;
Name = '';

% csv-Zeile in Bit-Array umwandeln
BA = serroblib_csvline2bits(csvline);

%% Bit-Array in Liste suchen
N = length(BA); % Anzahl Gelenke

% Daten laden
repopath=fileparts(which('serroblib_path_init.m'));
mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof');

% Filter für Typ des Roboters und DH-Parameter
num = 0;
for i = 1:size(l.BitArrays_Ndof, 1)
  % Zähle, wie viele Roboter dieser Gelenkreihenfolge existieren
  % (das erste Bit kennzeichnet den Gelenktyp)
  num = num + all(bitand(BA, 1) == bitand(l.BitArrays_Ndof(i,:), 1));
  
  % Prüfe ob alle Bits übereinstimmen (alle MDH-Parameter)
  if(all(BA == l.BitArrays_Ndof(i,:)))
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