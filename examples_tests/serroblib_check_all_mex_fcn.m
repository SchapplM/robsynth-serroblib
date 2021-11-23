% Rufe alle Matlab-Funktionen einmal auf und prüfe damit die Korrektheit
% der mex-Funktionen. Automatische Neu-Erzeugung, falls Syntax-Fehler.
% Dieses Skript ist dann nützlich, wenn die Schnittstellen im Robotik-Repo
% geändert wurden und noch mex-Dateien vorliegen, die mit der alten Version
% kompiliert wurden. Es werden auch die zugrunde liegenden m-Dateien geprüft.

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
    % Debug: Einzelnen Roboter prüfen.
%     if ~contains(Name, 'S6RRRRRR10V2'), continue; end
    isvar = l.AdditionalInfo(j,2); % Marker, ob Modellvariante
    Name_GenMdl = l.Names_Ndof{l.AdditionalInfo(j,3)};
    % Suche nach den m- und mex-Dateien
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
    mfilelist = dir(fullfile(tpl_dir, '*.m'));
    filelist = [mfilelist;mexfilelist; ];
    if isempty(filelist)
      % Nichts zu prüfen. Es gibt keine Dateien
      continue
    end
    fprintf('Prüfe Roboter %d/%d (%s) (%d mex-Dateien und %d m-Dateien liegen in tpl-Ordner %s)\n', ...
      j, length(l.Names_Ndof), Name, length(mexfilelist), length(mfilelist), tpl_dir)
    % Initialisiere Matlab-Klasse und setze auf Nutzung von M-Funktionen
    RS = serroblib_create_robot_class(Name);
    RS.gen_testsettings(true, true); % Zufallswerte für Parameter. Sonst kommen NaN-Fehler und überdecken die Syntax-Fehler
    % Gehe alle m- und mex-Dateien durch, die da sind. Fange mit m an.
    % Dadurch werden die Vorlagen-Funktionen meistens schon neu generiert.
    RS.fill_fcn_handles(false, false); RS_mex_status = false;
    for kk = 1:length(filelist)
      if contains(filelist(kk).name, '_mex')
        % Prüfe ab jetzt die mex-Dateien. Die m-Dateien sind fertig.
        if RS_mex_status == false
          RS.fill_fcn_handles(true, false);
          RS_mex_status = true;
        end
      end
      for retryiter = 1:3 % mehrfache Neuversuche zur Fehlerkorrektur
        recompile = false;
        % Einzelne Fälle für die mex-Dateien durchgehen und jeweils
        % Dummy-Aufruf der Funktionen, um Syntax-Fehler aufzudecken.
        if contains(filelist(kk).name, 'invkin_eulangresidual')
          try
            % Führe die Funktion aus
            % Gebe mehr als einen Startwert vor (neue Schnittstelle seit 2021-06)
            [~,~,~,Stats] = RS.invkin2(eye(3,4), rand(RS.NJ,2));
            % Prüfe, ob neue Ausgabe (seit 2021-06) da ist.
            tmp = Stats.coll;
            % Prüfe, ob Versionsausgabe existiert, und ob der Wert passt
            if Stats.version < 1 % hier wird die aktuelle Version eingetragen
              error('Version der Datei ist zu alt (%d).', Stats.version);
            end
            % Prüfe, ob Korrektur von Fehler bei Kollisionsprüfung da ist
            % Behoben ca. 2021-07; max/min mit Eingabe variabler Länge
            s = struct('avoid_collision_finish', true);
            [~,~,~,Stats] = RS.invkin2(eye(3,4), rand(RS.NJ,3), s);
          catch err
            if ~strcmp(err.identifier, 'MATLAB:svd:matrixWithNaNInf')
              recompile = true;
            end
          end
        end
        if contains(filelist(kk).name, 'invkin_traj')
          try
            % Führe die Funktion aus
           [~, ~, ~, ~, ~, Stats] = ...
             RS.invkin2_traj(zeros(2,6), zeros(2,6), zeros(2,6), [0;1], zeros(RS.NQJ,1));
            % Prüfe, ob Versionsausgabe existiert, und ob der Wert passt
            if Stats.version < 1 % hier wird die aktuelle Version eingetragen
              error('Version der Datei ist zu alt (%d).', Stats.version);
            end
          catch err
            if ~strcmp(err.identifier, 'MATLAB:svd:matrixWithNaNInf')
              recompile = true;
            end
          end
        end
        % Falls ein Fehler vorliegt, wird neu kompiliert oder generiert.
        if recompile
          warning('Fehler beim Aufruf von Funktion %s. Fehler: %s', ...
            filelist(kk).name, err.message);
          if retryiter == 2
            warning('Zweiter Versuch ohne Erfolg. Voraussichtlich sind die Vorlagen-Funktionen veraltet');
            serroblib_create_template_functions({Name}, false, false);
          elseif retryiter == 3
            error('Auch beim dritten Versuch kein Erfolg');
          end
          if contains(filelist(kk).name, '_mex')
            % Mex-Dateien werden neu kompiliert. m-Dateien nicht.
            [~,mexbasename] = fileparts(filelist(kk).name);
            matlabfcn2mex({mexbasename(1:end-4)}); % ohne "_mex"
          end
        else
          % Alles funktioniert. Keine Neuversuche notwendig.
          if retryiter > 1
            fprintf('Datei %s wurde korrigiert\n', filelist(kk).name);
          end
          break
        end % recompile
      end % retryiter
    end % kk (Mex-Dateien)
  end % j (Roboter-Modelle)
end % N (Gelenk-Anzahl)

