% Kompiliere alle (benötigten) Matlab-Funktionen für alle seriellen Robotermodelle
% 
% Nutzen dieses Skripts:
% * Nach Aufruf sind alle mex-Funktionen vorhanden und auf dem neusten
% Stand. Weniger Fehlermöglichkeiten durch veraltete mex-Dateien und
% schnellerer Aufruf anderer Skripte, da mex-Dateien nicht kompiliert
% werden müssen

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-03
% (C) Institut für Mechatronische Systeme, Universität Hannover
clear
clc
repopath=fileparts(which('serroblib_path_init.m'));

%% Funktionen für alle Modelle kompilieren
for N = 1:7
  % Liste zusammenstellen
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof');
  % alle Modelle generieren
  for iFK = 1:length(l.Names_Ndof)
    Name = l.Names_Ndof{iFK};
    fprintf('Modell %d/%d für %d FG: %s\n', ...
      iFK, length(l.Names_Ndof), N, Name);

    RS = serroblib_create_robot_class(Name);
    RS.fill_fcn_handles(true);
    RS.mex_dep(true); % Neu-Kompilierung erzwingen
  end
end
