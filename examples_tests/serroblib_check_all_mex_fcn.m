% Rufe alle Matlab-Funktionen einmal auf und prüfe damit die Korrektheit
% der mex-Funktionen. Automatische Neu-Erzeugung, falls Syntax-Fehler.
% Dieses Skript ist dann nützlich, wenn die Schnittstellen im Robotik-Repo
% geändert wurden und noch mex-Dateien vorliegen, die mit der alten Version
% kompiliert wurden.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
repopath=fileparts(which('serroblib_path_init.m'));
for N = 1:7 % Alle Gelenk-FG durchgehen
  % Liste zusammenstellen
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo', 'BitArrays_Ndof');
  fprintf('Prüfe %d Roboter mit %d Gelenken\n', length(l.Names_Ndof), N);
  for j = 1:length(l.Names_Ndof) % gehe alle Modelle durch
    Name = l.Names_Ndof{j};
    isvar = l.AdditionalInfo(j,2); % Marker, ob Modellvariante
    Name_GenMdl = l.Names_Ndof{l.AdditionalInfo(j,3)};
    % Suche nach den Mex-Dateien
    if ~isvar
      tpl_dir = fullfile(repopath, sprintf('mdl_%ddof', N), Name, 'tpl');
    else % Bei Varianten andere Ordnerstruktur
      tpl_dir = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, ...
        sprintf('tpl_%s', Name(length(Name_GenMdl)+1:end)));
      if ~isfolder(tpl_dir)
        tpl_dir = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, 'tpl');
      end
    end
    mexfilelist = dir(fullfile(tpl_dir, '*_mex.mex*'));
    if isempty(mexfilelist)
      % Nichts zu prüfen. Es gibt keine mex-Dateien
      continue
    end
    fprintf('Prüfe Roboter %d/%d (%s) (%d mex-Dateien liegen in tpl-Ordner %s)\n', ...
      j, length(l.Names_Ndof), Name, length(mexfilelist), tpl_dir)
    RS = serroblib_create_robot_class(Name);
    RS.fill_fcn_handles(true, false);
    % Gehe alle mex-Dateien durch, die da sind.
    for kk = 1:length(mexfilelist)
      for retryiter = 1:3 % mehrfache Neuversuche zur Fehlerkorrektur
        recompile = false;
        % Einzelne Fälle für die mex-Dateien durchgehen und jeweils
        % Dummy-Aufruf der Funktionen, um Syntax-Fehler aufzudecken.
        if contains(mexfilelist(kk).name, 'invkin_eulangresidual_mex')
          try
            RS.invkin2(NaN(3,4), zeros(RS.NJ,1));
          catch err
            recompile = true;
          end
        end
        if contains(mexfilelist(kk).name, 'invkin_traj_mex')
          try
            RS.invkin2_traj(NaN(2,6), NaN(2,6), NaN(2,6), [0;1], zeros(RS.NQJ,1));
          catch err
            recompile = true;
          end
        end
        % Falls ein Fehler vorliegt, wird neu kompiliert oder generiert.
        if recompile
          warning('Fehler beim Aufruf von Funktion %s. Fehler: %s', ...
            mexfilelist(kk).name, err.message);
          if retryiter == 2
            warning('Zweiter Versuch ohne Erfolg. Voraussichtlich sind die Vorlagen-Funktionen veraltet');
            serroblib_create_template_functions({Name}, false, false);
          elseif retryiter == 3
            error('Auch beim dritten Versuch kein Erfolg');
          end
          [~,mexbasename] = fileparts(mexfilelist(kk).name);
          matlabfcn2mex({mexbasename(1:end-4)}); % ohne "_mex"
        else
          % Alles funktioniert. Keine Neuversuche notwendig.
          if retryiter > 1
            fprintf('Mex-Datei %s wurde korrigiert\n', mexfilelist(kk).name);
          end
          break
        end % recompile
      end % retryiter
    end % kk (Mex-Dateien)
  end % j (Roboter-Modelle)
end % N (Gelenk-Anzahl)

