% Entferne Dateien aus dem Repo von Robotern, deren DH-Parameter
% aktualisiert wurden. Notwendig, da über Git-Befehle nicht alle Dateien
% aktualisiert werden.
% 
% Vorher ausführen: Skript was die Liste der Roboter erstellt. Z.B.:
% * remove_base_alignment_from_MDH.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

roblibpath=fileparts(which('serroblib_path_init.m'));

% Lade Liste mit Beinketten, die aktualisiert wurden. In Projektablage gespeichert.
papath = fileparts(which('robsynth_projektablage_path.m'));
listfile = fullfile(papath, '03_Entwicklung', 'Struktursynthese', ...
  '20220121_Aktualisierung_Beinketten_Basisausrichtung', ...
  'base_alignment_changed_list20220120_091120.txt');
LegNames_updated = readlines(listfile);
l = load(fullfile(roblibpath, 'serrob_list.mat'));

% Gehe alle seriellen Ketten aus der Liste durch und entferne jeweils den tpl-Ordner
for i = 1:length(LegNames_updated)
  Name_i = LegNames_updated{i};
  II = strcmp(l.Names, Name_i);
  if ~any(II), continue; end
  N = str2double(Name_i(2)); % Anzahl FG
  variantof = l.AdditionalInfo(II,3);
  isvariant = l.AdditionalInfo(II,2);
  if isvariant
    % Prüfe, ob der Ordner des Variantenmodells existiert
    Name_GenMdl = l.Names{variantof};
    if ~contains(Name_GenMdl, Name_i), error('Name passt nicht'); end
    fcn_dir = fullfile(roblibpath, sprintf('mdl_%ddof', N), Name_GenMdl, ...
      sprintf('tpl_%s', Name_i(length(Name_GenMdl)+1:end)));
  else
    fcn_dir = fullfile(roblibpath, sprintf('mdl_%ddof', N), Name_i, 'tpl');
  end
  if exist(fcn_dir, 'file')
    fprintf('Verzeichnis %s gelöscht. Gehört zu aktualisierter Kinematik.\n', fcn_dir);
    rmdir(fcn_dir, 's');
  end
end
