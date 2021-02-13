% Füge die Modelle gegebener Roboter zum Pfad hinzu
% 
% Eingabe:
% Names
%   Cell-Array mit Namen der Robotermodelle, deren Funktionen zum
%   Matlab-Pfad hinzugefügt werden sollen.
% 
% Ausgabe:
% added
%   Bool-Array. "true", falls für den Eintrag aus `Names` der Pfad geändert
%   wurde. "false", falls der Roboter schon vollständig im Pfad war.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Leibniz Universität Hannover

function added = serroblib_addtopath(Names)
added = false(length(Names),1);
if ~iscell(Names)
  error('serroblib_addtopath: Eingegebene Roboternamen müssen ein Cell-Array sein');
end
repopath=fileparts(which('serroblib_path_init.m'));
for i = 1:length(Names)
  Name = Names{i};
  N = str2double(Name(2)); % Nehme Namensschema SxRRRR mit x=Anzahl FG

  % Prüfe, ob Code für dieses Modell existiert
  fcn_dir = fullfile(repopath, sprintf('mdl_%ddof', N), Name, 'hd');
  if exist(fcn_dir, 'file')
    added(i) = added(i) | addpath2(fcn_dir);
    % Falls Ordner für Vorlagen-Funktionen nicht existiert, erstelle neu
    tpl_dir = fullfile(repopath, sprintf('mdl_%ddof', N), Name, 'tpl');
    if ~exist(tpl_dir, 'file')
      serroblib_create_template_functions({Name});
      % Pfad kann durch create_template_functions wieder gelöscht werden.
      % Laufzeitverhalten durch mehrfache Rekursion nicht genau bestimmbar.
      % Daher Zeile behalten (diese Funktion kann von template-Funktion
      % aufgerufen werden.
      added(i) = added(i) | addpath2(fcn_dir);
    end
    added(i) = added(i) | addpath2(tpl_dir);
  else
    % Prüfe, ob Modell eine Variante ist, für die der Code des Hauptmodells
    % existiert
    mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
    l = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo');
    isvariant = l.AdditionalInfo(strcmp(l.Names_Ndof,Name),2);
    variantof = l.AdditionalInfo(strcmp(l.Names_Ndof,Name),3);
    if isvariant
      % Prüfe, ob der Ordner des Variantenmodells existiert
      Name_GenMdl = l.Names_Ndof{variantof};
      fcn_dir_var = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, ...
        sprintf('hd_%s', Name(length(Name_GenMdl)+1:end)));
      if exist(fullfile(fcn_dir_var, sprintf('%s_structural_kinematic_parameters.m', Name)), 'file')
        % Es wurde schon Code spezifisch für diese Variante generiert
        added(i) = added(i) | addpath2(fcn_dir_var);
        % Vorlagen-Funktionen sollten dann im selben Format vorliegen
        tpl_dir_var = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, ...
          sprintf('tpl_%s', Name(length(Name_GenMdl)+1:end)));
        if ~exist(tpl_dir_var, 'file')
          serroblib_create_template_functions({Name});
          added(i) = added(i) | addpath2(fcn_dir_var); % Pfad wird durch create_template_functions wieder gelöscht
        end
        added(i) = added(i) | addpath2(tpl_dir_var);
      else % Kein spezifischer Code da
        % Prüfe, ob Funktionen zur Umwandlung von Allgemein zu Variante da
        % ist
        var_dir = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, 'var');
        pkinconvfile = sprintf('%s_pkin_gen2var.m', Name);
        if ~exist(fullfile(var_dir,pkinconvfile), 'file')
          warning('Varianten-Funktionen für %s in %s existieren nicht. Spezifischer Code in %s existiert auch nicht.', ...
            Name_GenMdl, var_dir, fcn_dir_var);
        else
          % Varianten-Funktion existiert. Füge Pfad hinzu
          added(i) = added(i) | addpath2(var_dir);
        end
        fcn_dir_gen = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, 'hd');
        if exist(fcn_dir_gen, 'file')
          added(i) = added(i) | addpath2(fcn_dir_gen);
        else
          warning('Verzeichnis für allgemeinen Code %s existiert nicht. Es wird aber darauf verwiesen.', fcn_dir_gen);
        end
        tpl_dir_gen = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, 'tpl');
        if ~exist(tpl_dir_gen, 'file')
          serroblib_create_template_functions({Name_GenMdl});
          % Pfad kann durch create_template_functions wieder gelöscht werden.
          added(i) = added(i) | addpath2(fcn_dir_gen);
        end
        added(i) = added(i) | addpath2(tpl_dir_gen);
      end
    else
      % Ist keine Variante und Verzeichnis nicht da.
      warning('Code-Verzeichnis %s für %s existiert nicht', fcn_dir, Name);
    end
  end
end
end

function added = addpath2(new_dir)
% Füge einen Ordner zum Pfad hinzu, wenn dieser bislang noch gefehlt hat.
% Daduch wird die Reihenfolge eines bereits im Pfad befindlichen Ordners
% nicht mehr geändert.
if ~contains(path(), [new_dir, pathsep()])
  addpath(new_dir);
  added = true;
else
  added = false;
end
end