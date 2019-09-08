% Entferne die Modelle gegebener Roboter vom Matlab-Pfad

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-10
% (C) Institut für mechatronische Systeme, Universität Hannover

function serroblib_removefrompath(Names)
repopath=fileparts(which('serroblib_path_init.m'));
exp_var = '^S(\d)([RP]+)(\d+)V(\d+)$'; % Format "S3RRR1V1" für Varianten

for i = 1:length(Names)
  n = Names{i};
  N = str2double(n(2)); % Nehme Namensschema SxRRRR mit x=Anzahl FG
  % Prüfe, ob Modellname eine Variante ist
  [tokens_var, ~] = regexp(n,exp_var,'tokens','match');
  if ~isempty(tokens_var) % Variante
    % Variante: Entferne den speziellen Variantenordner aus dem Pfad
    Name_GenMdl = sprintf('S%s%s%s', tokens_var{1}{1}, tokens_var{1}{2}, tokens_var{1}{3});
    fcn_dir_var = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, ...
      sprintf('hd_V%s', tokens_var{1}{4}));
    tpl_dir_var = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, ...
      sprintf('tpl_V%s', tokens_var{1}{4}));
    if exist(fcn_dir_var, 'file')
      rmpath(fcn_dir_var);
    end
    if exist(tpl_dir_var, 'file')
      rmpath(tpl_dir_var);
    end
    % Entferne den Ordner des Hauptmodells
    fcn_dir_gen = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, 'hd');
    tpl_dir_gen = fullfile(repopath, sprintf('mdl_%ddof', N), Name_GenMdl, 'tpl');
    if ~exist(fcn_dir_var, 'file') && exist(fcn_dir_gen, 'file')
      rmpath(fcn_dir_gen);
    end
    if ~exist(tpl_dir_var, 'file') && exist(tpl_dir_gen, 'file')
      rmpath(tpl_dir_gen);
    end
  else % Keine Variante
    % Entferne den Dateien-Ordner dieses expliziten (allgemeinen
    % Robotermodells aus dem Pfad
    fcn_dir = fullfile(repopath, sprintf('mdl_%ddof', N), n, 'hd');
    tpl_dir = fullfile(repopath, sprintf('mdl_%ddof', N), n, 'tpl');
    if exist(fcn_dir, 'file')
      rmpath(fcn_dir);
    end
    if exist(tpl_dir, 'file')
      rmpath(tpl_dir);
    end
  end
end