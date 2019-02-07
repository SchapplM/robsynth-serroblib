% Finde Robotermodelle in der Datenbank mit gegebenen Filtern
% 
% Eingabe:
% N
%   Anzahl der Freiheitsgrade des Roboters
% EE_FG [1x9] oder [1x6]
%   Vektor mit 1/0 für Belegung ob EE-FG aktiv ist
% EE_FG_Mask [1x9] oder [1x6]
%   Maske, die festlegt ob die FG exakt wie in `EE_FG` sind, oder ob auch
%   gesperrte FG wirklich nicht beweglich sind
% 
% Rückgabe:
% Indizes
%   Vektor mit Zähl-Indizes für zur Auswahl der passenden Roboter in der
%   Liste von allen Robotern mit `N` FG (Datei S3list.mat für 3FG).

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function [Indizes] = serroblib_filter_robots(N, EE_FG, EE_FG_Mask)

if length(EE_FG) == 6 && length(EE_FG_Mask) == 6
  % Altes Format: Ignoriere die Erweiterung um die Euler-Winkel des EE
  % Durch die Maske sind die Euler-Winkel egal
  EE_FG_Mask = [EE_FG_Mask, false(1,3)];
  % Setze gewünschte Euler-Winkel auf Winkelgeschwindigkeit (passt bei 1R
  % und 3R; nicht bei 2R). Aber da muss man sowieso die Euler-Winkel nehmen
  EE_FG = [EE_FG, EE_FG(4:6)];
end

% Umwandlung in Cell-Array-Format
EE_FG_cell = cell(1,9);
for i = 1:9
  EE_FG_cell{i} = sprintf('%d', EE_FG(i));
end
EE_FG_BA = serroblib_csvline2bits_EE(EE_FG_cell);

% Bitmaske als Variable speichern
EE_FG_Mask_bin = uint16(Inf);
for i = 1:9
  if EE_FG_Mask(i)
    EE_FG_Mask_bin = bitset(EE_FG_Mask_bin, i);
  end
end

repopath=fileparts(which('serroblib_path_init.m'));

Indizes = [];
b = 1;
mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_EEdof0');

% Bit-Maske für EE-FG
% Suche alle Roboter für FG N heraus
for i = 1:size(l.BitArrays_Ndof, 1)
  if (all( bitand(EE_FG_BA,EE_FG_Mask_bin) == bitand(l.BitArrays_EEdof0(i,:),EE_FG_Mask_bin) ))
    % Roboter mit passenden EE-FG
    Indizes(b) = i; %#ok<AGROW>
    b = b+1;      
  end
end