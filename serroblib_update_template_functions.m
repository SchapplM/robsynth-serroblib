% Prüfe die Template-Funktionen für einen bestimmten Roboter.
% Aktualisiere die Funktionen, falls sie veraltet sind.
% 
% Eingabe:
% Names {1 x n} cell array
%   Namen der Roboter, für die die Funktionen aktualisiert werden.
%   Falls leer gelassen: Alle.
% verbosity
%   Grad der Textausgabe (0=aus, 1=Fortschritt)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function serroblib_update_template_functions(Names, verbosity)
orig_state = warning('off', 'all'); % Warnungen temporär unterdrücken
%% Prüfe Eingabe
repopath=fileparts(which('serroblib_path_init.m'));
% Stelle Liste aller Roboter zusammen
mdllistfile_alldof = fullfile(repopath, 'serrob_list.mat');
l = load(mdllistfile_alldof, 'Names', 'AdditionalInfo');
if nargin < 1 || isempty(Names) % keine Eingabe einer Liste. Nehme alle.
  Names = l.Names;
end
if nargin < 2
  verbosity = 0;
end
%% Bestimme aktuelle Version der jeweiligen Vorlagen-Funktion
filelist = {'robot_invkin_eulangresidual.m.template', 'robot_invkin_traj.m.template'};
fileversions = struct('robot_invkin_eulangresidual', 0, 'robot_invkin_traj', 0);
kinematics_dir = fullfile(fileparts(which('robotics_toolbox_path_init.m')), ...
  'kinematics');
for i = 1:length(filelist)
  [tokens,~] = regexp(fileread(fullfile(kinematics_dir, filelist{i})), ...
    '''version'', (\d+)', 'tokens', 'match');
  key = filelist{i}(1:end-11); % Entferne Endung .m.template
  fileversions.(key) = str2double(tokens{1}{1});
end

%% Gehe alle Dateien durch und prüfe die Version
for j = 1:length(Names)
  Name = Names{j};
  AdditionalInfo_j = l.AdditionalInfo(strcmp(Name, l.Names),:);
  if isempty(AdditionalInfo_j)
    warning('Roboter %s nicht in Datenbank gefunden', Name);
  end
  isvar = AdditionalInfo_j(2); % Marker, ob Modellvariante
  Name_GenMdl = l.Names{AdditionalInfo_j(3)};
  N = str2double(Name(2));
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
  if verbosity
    fprintf('Prüfe Roboter %d/%d (%s) (%d mex-Dateien und %d m-Dateien liegen in tpl-Ordner %s)\n', ...
      j, length(Names), Name, length(mexfilelist), length(mfilelist), tpl_dir)
  end
  % Initialisiere Matlab-Klasse und setze auf Nutzung von M-Funktionen
  RS = serroblib_create_robot_class(Name);
  RS.gen_testsettings(true, true); % Zufallswerte für Parameter. Sonst kommen NaN-Fehler und überdecken die Syntax-Fehler
  RS.qlim = repmat([-1,1],RS.NQJ,1);
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
          s = struct('n_max', 1, 'retry_limit', -1); % keine Ausführung
          [~,~,~,Stats] = RS.invkin2(eye(3,4), rand(RS.NJ,2), s);
          % Prüfe, ob neue Ausgabe (seit 2021-06) da ist.
          if Stats.version < fileversions.robot_invkin_eulangresidual
            error('Version der Datei ist zu alt (%d). Aktuell: %d', ...
              Stats.version, fileversions.robot_invkin_eulangresidual);
          end
        catch err
          recompile = true;
        end
      end
      if contains(filelist(kk).name, 'invkin_traj')
        try
          % Führe die Funktion aus
         [~, ~, ~, ~, ~, Stats] = ...
           RS.invkin2_traj(zeros(2,6), zeros(2,6), zeros(2,6), [0;1], zeros(RS.NQJ,1));
          % Prüfe, ob Versionsausgabe existiert, und ob der Wert passt
          if Stats.version < fileversions.robot_invkin_traj
            error('Version der Datei ist zu alt (%d). Aktuell: %d', ...
              Stats.version, fileversions.robot_invkin_traj);
          end
        catch err
          if ~strcmp(err.identifier, 'MATLAB:svd:matrixWithNaNInf')
            recompile = true;
          end
        end
      end
      % Falls ein Fehler vorliegt, wird neu kompiliert oder generiert.
      if recompile
        if verbosity
          fprintf('Fehler beim Aufruf von Funktion %s. Fehler: %s\n', ...
            filelist(kk).name, err.message);
        end
        if retryiter == 2
          if verbosity
            fprintf(['Zweiter Versuch ohne Erfolg. Voraussichtlich sind ', ...
              'die Vorlagen-Funktionen veraltet\n']);
          end
          serroblib_create_template_functions({Name}, false, false);
        elseif retryiter == 3
          warning(orig_state);
          error('Auch beim dritten Versuch kein Erfolg');
        end
        if contains(filelist(kk).name, '_mex')
          % Mex-Dateien werden neu kompiliert. m-Dateien nicht.
          [~,mexbasename] = fileparts(filelist(kk).name);
          matlabfcn2mex({mexbasename(1:end-4)}); % ohne "_mex"
        end
      else
        % Alles funktioniert. Keine Neuversuche notwendig.
        if retryiter > 1 && verbosity
          fprintf('Datei %s wurde korrigiert\n', filelist(kk).name);
        end
        break
      end % recompile
    end % retryiter
  end % kk (Mex-Dateien)
end % j (Roboter-Modelle)
warning(orig_state);