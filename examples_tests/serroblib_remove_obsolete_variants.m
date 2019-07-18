% Entferne Varianten von Modellen, die nicht mehr existieren

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-06
% (C) Institut für Mechatronische Systeme, Universität Hannover

clc
clear

serroblib_gen_bitarrays;

%% Hauptmodelle entfernen
% Beispielhafte Liste von Modellen (Lösch-Liste vom 14.06.2019)
Names_niO = {...
    'S6RPRRRP7' , ...
    'S6RPRRRR1' , ...
    'S6RPRRRR6' , ...
    'S6RPRRRR7' , ...
    'S6RRPRRP4' , ...
    'S6RRPRRP12', ...
    'S6RRPRRR3' , ...
    'S6RRPRRR6' , ...
    'S6RRPRRR11', ...
    'S6RRRPRP9' , ...
    'S6RRRPRR3' , ...
    'S6RRRPRR10', ...
    'S6RRRRPP6' , ...
    'S6RRRRPR2' , ...
    'S6RRRRPR3' , ...
    'S6RRRRPR8' , ...
    'S6RRRRRP1' , ...
    'S6RRRRRP6' , ...
    'S6RRRRRR1' , ...
    'S6RRRRRR4'};
  
for i = 1:length(Names_niO)
  success = serroblib_remove_robot(Names_niO{i});
  if success
    fprintf('Roboter %s erfolgreich gelöscht\n', Names_niO{i});
  else
    fprintf('Roboter %s konnte nicht gelöscht werden\n', Names_niO{i});
  end
end

%% Varianten entfernen
for N = 1:7
  serroblibpath=fileparts(which('serroblib_path_init.m'));
  mdllistfile_Ndof = fullfile(serroblibpath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l_serrob = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo');
  
  II_obs =( l_serrob.AdditionalInfo(:,2) == 1) & (l_serrob.AdditionalInfo(:,3) == 0);
  for i = find(II_obs)'
    Name_i = l_serrob.Names_Ndof{i};
    success = serroblib_remove_robot(Name_i);
    fprintf('Variante %s ohne Hauptmodell aus Datenbank entfernt\n', Name_i);
  end
end