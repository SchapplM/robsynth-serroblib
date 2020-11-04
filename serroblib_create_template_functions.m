% Erstelle die Vorlagen-Funktionen aller Robotermodelle in der Bibliothek
% 
% Eingabe:
% Names
%   Cell array mit Namen der zu erstellenden Robotermodelle
%   Optional: Wenn nicht angegeben, werden alle erstellt.
% skip_existing
%   true: Bereits existierende tpl-Ordner überspringen (Standard)
%   false: Alle Dateien immer neu erstellen
% mex_results
%   true: Die aus Vorlagen generierten Funktionen werden zum testen direkt
%   kompiliert.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07
% (C) Institut für Mechatronische Systeme, Universität Hannover

function serroblib_create_template_functions(Names, skip_existing, mex_results)

%% Initialisierung
old_dir = pwd();
% Prüfe Eingabeargument
repopath=fileparts(which('serroblib_path_init.m'));
if nargin < 1 || isempty(Names)
  % Stelle Liste aller Roboter zusammen
  Names = {};
  for N = 1:7
    mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
    if ~exist(mdllistfile_Ndof, 'file')
      continue
    end
    % Lade in .mat gespeicherte Datenbank
    l = load(mdllistfile_Ndof, 'Names_Ndof');
    Names = {Names{:}, l.Names_Ndof{:}}; %#ok<CCAT>
  end
end
if nargin < 2
  skip_existing = true;
end
if nargin < 3
  mex_results = false;
end
% Kopier-/Ersetz-Befehl vorbereiten (zwei Alternativen
if ispc()
  % Der sed-Befehl unter Windows ist möglich (wenn Git for Windows
  % installiert ist), ist aber sehr langsam und wird nicht benutzt.
  sedcmd = '"C:\\Program Files\\Git\\usr\\bin\\sed.exe" -i ''s/%s/%s/g'' %s';
  copycommand = 'linewise';
else
  sedcmd = 'sed -i ''s/%s/%s/g'' %s';
  copycommand = 'sed';
end
%% Alle Vorlagen-Funktionen aus HybrDyn-Repo kopieren
% Alle Funktionen dieser Liste werden roboterspezifisch aus der Liste
% erstellt. Die Anpassung sind nur geringfügig und ermöglichen Kompilierung
function_list_copy_hybrdyn = {...
  'gravload_floatb_eulxyz_nnew_vp1.m', ...
  'jacobig_mdh_num.m', ...
  'jacobigD_mdh_num.m', ...
  'jacobig_cutforce_mdh_num.m', ...
  'inertiaJ_nCRB_vp1.m', ...
  'invdyn_floatb_eulxyz_nnew_vp1.m', ...
  'invdyn_floatb_eulxyz_nnew_vp2.m', ...
  'gravload_floatb_eulxyz_nnew_vp1.m', ...
  'invdynJ_fixb_mdp_slag_vp_traj.m', ...
  'invdynJ_fixb_mdp_slag_vr_traj.m', ...
  'invdynJ_fixb_regmin_slag_vp_traj.m'};
function_list_copy_robotics = {...
  {'kinematics', 'constr2.m'}, ...
  {'kinematics', 'constr2grad.m'}, ...
  {'kinematics', 'invkin_traj.m'}, ...
  {'kinematics', 'invkin_eulangresidual.m'}};

% Pfad zur Maple-Dynamik-Toolbox, in der einige Vorlagen liegen
mrp = fileparts(which('hybrdyn_path_init.m'));
if isempty(mrp)
  warning('Die HybridDyn-Toolbox muss im Pfad sein (siehe README.MD)');
  return
end
% Alle Vorlagen-Funktionen aus Maple-Repo in Vorlagen-Ordner kopieren
for tmp = function_list_copy_hybrdyn
  tplf = tmp{1};
  copyfile(fullfile(mrp, 'robot_codegen_scripts', 'templates_num', ['robot_',tplf,'.template']), ...
           fullfile(repopath, 'template_functions') );
end
rtp = fileparts(which('robotics_toolbox_path_init.m'));
if isempty(rtp)
  warning('Die Robotik-Toolbox muss im Pfad sein (siehe README.MD)');
  return
end
for tmp = function_list_copy_robotics
  tplf = tmp{1};
  copyfile(fullfile(rtp, tplf{1}, ['robot_',tplf{2},'.template']), ...
           fullfile(repopath, 'template_functions') );
end

% Generiere die Liste der Vorlagen-Funktionen aus den tatsächlich
% existierenden Dateien. Dadurch können durch Benutzer hinzugefügte
% Funktionen auch generiert werden
fl_tmp = dir(fullfile(repopath, 'template_functions', '*.template'));
function_list = cell(1,length(fl_tmp));
for i = 1:length(fl_tmp)
  % Präfix "robot_" und Suffix ".template" entfernen
  function_list{i} = strrep(fl_tmp(i).name(7:end), '.template', '');
end
%% Gehe alle Modellnamen durch und erstelle alle Vorlagen-Funktionen
for i = 1:length(Names)
  Name_i = Names{i};
  N = str2double(Name_i(2)); % Anzahl FG
  % Daten des Robotermodells laden
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo');
  I_robot = find(strcmp(l.Names_Ndof,Name_i));
  if isempty(I_robot)
    error('Modell %s ist nicht in der Datenbank', Name_i);
  end
  addinfo = l.AdditionalInfo(I_robot,:);
  isvariant = addinfo(2);
  variantof = addinfo(3);
  hascode   = addinfo(4);
	% Prüfe, ob Roboter einen Code-Ordner hat. Wenn nicht, ist die Erstellung
	% nicht sinnvoll (weil die robot_env.sh nicht existiert)
  if hascode ~= 1
    if isvariant
      % Stattdessen die Funktionen für das allgemeine Modell erzeugen
      Name_Gen_i = l.Names_Ndof{variantof};
      serroblib_create_template_functions({Name_Gen_i}, skip_existing, mex_results);
    end
    continue
  end
  % Bestimme Ziel-Ordner des Codes für die Vorlagen-Funktionen
  if isvariant
    % Prüfe, ob der Ordner des Variantenmodells existiert
    Name_GenMdl = l.Names_Ndof{variantof};
    fcn_dir = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, ...
      sprintf('tpl_%s', Name_i(length(Name_GenMdl)+1:end)));
    def_file = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, ...
      sprintf('hd_%s', Name_i(length(Name_GenMdl)+1:end)), sprintf('robot_env_%s.sh', Name_i));
  else
    fcn_dir = fullfile(repopath, sprintf('mdl_%ddof', N), Name_i, 'tpl');
    def_file = fullfile(repopath, sprintf('mdl_%ddof', N), Name_i, 'hd', sprintf('robot_env_%s.sh', Name_i));
  end
  
  % Definitionsdatei lesen
  filetext = fileread(def_file);
  % Platzhalter-Ausdrücke für diesen Roboter erhalten (aus robot_env.sh)
  subsexp_array = {'robot_name', 'RN', {}; ...
                   'robot_NQJ', 'NQJ', {}; ...
                   'robot_NJ',  'NJ',  {}; ...
                   'robot_NL',  'NL',  {}; ...
                   'robot_NMPVFIXB',  'NMPVFIXB',  {}; ...
                   'robot_NMPVFLOATB',  'NMPVFLOATB',  {}; ...
                   'robot_NKP',  'NKP',  {}; ...
                   'robot_KP',  'KPDEF',  {}; ...
                   'robot_NTAUJFIXBREGNN',  'NTAUJFIXBREGNN',  {}};
  for ii = 1:size(subsexp_array,1)
    expr = [subsexp_array{ii,1}, '=(.*)'];
    tokens = regexp(filetext,expr,'tokens','dotexceptnewline');
    subsexp_array(ii,3) = strrep(tokens{1},'"','');
    subsexp_array(ii,3) = strrep(subsexp_array(ii,3),sprintf('\r'),''); % CR-Zeichen entfernen (Windows/Linux-Problem)
    if strcmp(subsexp_array{ii,2}, 'KPDEF')
      subsexp_array{ii,3} = sprintf('pkin: %s', subsexp_array{ii,3});
    end
  end
  subsexp_array{10,2} = 'VERSIONINFO';
  subsexp_array{10,3} = 'Generated in SerRobLib from HybrDyn-Template';
  
  % Kopiere alle Vorlagen-Funktionen an die Ziel-Orte und Ersetze die
  % Platzhalter-Ausdrücke
  mkdirs(fcn_dir);  % Ordner existiert noch nicht. Neu erstellen.
  cd(fcn_dir); % In Ordner wechseln für kürzeren sed-Befehl (und zum Finden der Dateien)
  function_list_mex = {};
  num_files_written = 0;
  for tmp = function_list
    tplf = tmp{1};
    file1=fullfile(repopath, 'template_functions', ['robot_',tplf,'.template']);
    file2=fullfile(fcn_dir, [Name_i,'_',tplf]);
    function_list_mex = [function_list_mex(:); [Name_i,'_',tplf(1:end-2)]];
    if skip_existing && exist(file2, 'file')
      continue % Keine Datei erzeugen. Nur mex-Liste erstellen.
    end
    % Ersetzungsausdruck für Dateinamen vorbereiten
    subsexp_array{11,2} = 'FN';
    [~,subsexp_array{11,3},~] = fileparts(file2);
    % Kopieren und Ausdrücke ersetzen
    if strcmp(copycommand, 'sed')
      % Variante für Linux: Vorlagen-Datei kopieren und anschließend
      % Text-Ersetzung mit sed-Befehl durchführen
      copyfile(file1, file2);
      for ii = 1:size(subsexp_array,1)
        sedcmd_ii = sprintf(sedcmd, ['%',subsexp_array{ii,2},'%'], subsexp_array{ii,3}, [Name_i,'_',tplf]);
        system(sedcmd_ii);
      end
    elseif strcmp(copycommand, 'linewise')
      % Variante für Windows: Text ersetzen mit sed ist sehr langsam.
      % Nutze daher die Zeilenweise kopier-Funktion von Matlab.
      fid1 = fopen(file1);
      fid2 = fopen(file2, 'w');
      tline = fgetl(fid1);
      while ischar(tline)
        % Zeile weiterverarbeiten: Platzhalter-Ausdrücke ersetzen
        for ii = 1:size(subsexp_array,1)
          tline = strrep(tline, ['%',subsexp_array{ii,2},'%'], subsexp_array{ii,3});
        end
        fwrite(fid2, [tline, newline()]); % Zeile in Zieldatei schreiben
        tline = fgetl(fid1); % nächste Zeile
      end
      fclose(fid1);fclose(fid2);
    else
      error('Befehl nicht definiert');
    end
    % fprintf('%d/%d: Vorlagen-Funktion %s erstellt.\n', i, length(Names), tplf);
    
    % Prüfe, ob die Erstellung erfolgreich war (falls die Code-Generierung
    % nicht vollständig durchgeführt wurde, gehen nicht alle Funktionen)
    fid2 = fopen(file2, 'r');
    tline = fgetl(fid2);
    function_invalid = false;
    while ischar(tline)
      if contains(tline, 'NOTDEFINED') % Wird in HybrDyn eingesetzt, falls Ausdruck nicht gefunden.
        function_invalid = true;
        break;
      end
      tline = fgetl(fid2);
    end
    fclose(fid2);
    if function_invalid
      fprintf('%d/%d: Datei %s konnte nicht erzeugt werden (Einige Variablen nicht definiert)\n', ...
        i, length(Names), tplf);
      delete(file2);
      continue;
    end
    num_files_written = num_files_written + 1;
  end
  fprintf('%d/%d: Vorlagen-Funktionen für %s erstellt (%d/%d).\n', i, ...
    length(Names), Name_i, num_files_written, length(function_list));
  
  % Prüfe, ob mex-Datei im tpl-Ordner liegt. Wenn nicht, ist wahrscheinlich
  % noch eine alte Version vorhanden. Diese Funktion ist sinnvoll zur
  % Aktualisierung der Repo-Version auf das tpl-Format
  serroblib_addtopath({Name_i});
  for tmp = function_list
    [~,f_basename] = fileparts([Name_i, '_', tmp{1}]);
    [dir_mexfcn, ~, mex_ext] = fileparts(which(sprintf('%s_mex', f_basename)));
    if isempty(dir_mexfcn)
      % mex-Datei existiert nicht (in irgendeinem Ordner im Pfad).
      continue
    elseif ~strcmp(dir_mexfcn, fcn_dir)
      if exist(fullfile(fcn_dir,[f_basename, '_mex', mex_ext]), 'file')
        % Die Mex-Datei liegt im Matlab-Pfad nicht nur im vorgesehenen
        % tpl-Ordner, sondern zusätzlich woanders. Lösche die "falsche"
        % mex-Datei, die wahrscheinlich veraltet ist
        delete(fullfile(dir_mexfcn, [f_basename, '_mex', mex_ext]));
        fprintf('Doppelte mex-Funktion %s in % s gelöscht\n', f_basename, dir_mexfcn);
      else
        % Die mex-Datei existiert nur im "falschen" Ordner. Verschiebe in
        % tpl-Ordner
        movefile(fullfile(dir_mexfcn, [f_basename, '_mex', mex_ext]), ...
                 fullfile(fcn_dir,    [f_basename, '_mex', mex_ext]));
        fprintf('Mex-Datei %s_mex%s an richtigen Ort (%s) verschoben (von %s)\n', ...
          f_basename, mex_ext, fcn_dir, dir_mexfcn);
      end
    else
      % Die mex-Datei existiert im vorgesehenen Ordner. Alles i.O.
      continue
    end
  end
  serroblib_removefrompath({Name_i});
  % Testen: Kompilieren aller Funktionen im Zielordner
  if mex_results
    serroblib_addtopath({Name_i})
    cd(fcn_dir)
    matlabfcn2mex(function_list_mex);
    serroblib_removefrompath({Name_i})
  end
end
% Zurückwechseln in vorheriges Verzeichnis
cd(old_dir);
