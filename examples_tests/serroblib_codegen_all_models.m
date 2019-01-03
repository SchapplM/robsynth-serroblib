% Generieren des Matlab-Codes für alle seriellen Robotermodelle
% 
% Nutzen dieses Skripts:
% * Neu-Generierung aller Matlab-Funktionen, z.B. wenn Änderungen in der
%   HybrDyn-Toolbox erfolgt sind.
% * Generieren von Matlab-Funktionen für Roboter, die neu zur
%   Modelldatenbank hinzugefügt wurden

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

repopath=fileparts(which('serroblib_path_init.m'));
%% Alle Robotermodell initialisieren
serroblib_gen_bitarrays
serroblib_create_robot_csv_all

%% Code für alle Modelle generieren
for N = 1:7
  % Liste zusammenstellen
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof');
  % alle Modelle generieren
  for iFK = 1:length(l.Names_Ndof)
    Name = l.Names_Ndof{iFK};
    serroblib_generate_mapleinput({Name})
    serroblib_generate_code({Name}, true)
    % serroblib_addtopath({Name})
  end
end

%% Code für ein bestimmtes Modell neu generieren
% Beispiel: S6RRPRRR14 (z.B. UPS-Kette)
Name = 'S6RRPRRR14';
serroblib_generate_mapleinput({Name})
serroblib_generate_code({Name}, true)
