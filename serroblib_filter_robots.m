% Finde Robotermodelle in der Datenbank mit gegebenen Filtern
% 
% Eingabe:
% N
%   Anzahl der Freiheitsgrade des Roboters
% EE_FG [1x9] oder [1x6]
%   Vektor mit 1/0 für Belegung ob EE-FG aktiv ist
% EE_FG_Mask [1x9] oder [1x6]
%   Maske, die festlegt ob die FG exakt wie in `EE_FG` sind (Bit auf 1), 
%   Bit 1: Gesperrter FG in `EE_FG` (auf 0) muss auch gesperrt sein, FG in
%   `EE_FG` (auf 1) muss auch frei sein
%   Bit 0: Bit in `EE_FG` wird nicht beachtet (z.B. zur Auswahl aller
%   Beinketten, die 3T haben, egal ob 0R,1R,2R,3R)
% 
% Rückgabe:
% IndZ
%   Vektor mit Zähl-Indizes für zur Auswahl der passenden Roboter in der
%   Liste von allen Robotern mit `N` FG (Datei S3list.mat für 3FG).
% IndB
%   Gleicher Vektor als Binärindex (1, falls entsprechender Roboter passend)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [IndZ, IndB] = serroblib_filter_robots(N, EE_FG, EE_FG_Mask)

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
EE_FG_Mask_bin = uint16(0);
for i = 1:9
  if EE_FG_Mask(i)
    EE_FG_Mask_bin = bitset(EE_FG_Mask_bin, i);
  end
end

repopath=fileparts(which('serroblib_path_init.m'));

mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_EEdof0');
IndB = false(length(l.Names_Ndof),1);

% Bit-Maske für EE-FG
% Suche alle Roboter für FG N heraus
for i = 1:size(l.BitArrays_Ndof, 1)
  if (all( bitand(EE_FG_BA,EE_FG_Mask_bin) == bitand(l.BitArrays_EEdof0(i,:),EE_FG_Mask_bin) ))
    % Roboter mit passenden EE-FG
    IndB(i) = true;
  end
  % Debug:
%   fprintf('%s: Soll: %s. Ist: %s. Maske %s. Passend: %d\n', l.Names_Ndof{i}, fliplr(dec2bin(EE_FG_BA,9)), ...
%     fliplr(dec2bin(l.BitArrays_EEdof0(i,:),9)), fliplr(dec2bin(EE_FG_Mask_bin,9)), IndB(i));
end
IndZ = find(IndB);