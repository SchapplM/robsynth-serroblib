% Kompiliere alle Matlab-Funktionen aller Robotermodelle in der Bibliothek
% 
% Eingabe:
% Names
%   Cell array mit Namen der zu kompilierenden Robotermodelle
%   Optional: Wenn nicht angegeben, werden alle kompiliert.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für Mechatronische Systeme, Universität Hannover

function serroblib_compile_functions(Names)

% Prüfe Eingabeargument
repopath=fileparts(which('serroblib_path_init.m'));
if nargin == 0
  % Stelle Liste aller Roboter zusammen
  Names = {};
  for N = 1:7
    mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
    if ~exist(mdllistfile_Ndof, 'file')
      continue
    end
    % Lade in .mat gespeicherte Datenbank
    l = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo');
    % Nur Roboter zur Liste hinzufügen, bei denen ein Code-Ordner vorliegt.
    Names = {Names{:}, l.Names_Ndof{l.AdditionalInfo(:,4)~=0}}; %#ok<CCAT>
  end
end

% Gehe alle Modellnamen durch und kompiliere alle m-Dateien in dem
% jeweiligen Ordner.
for i = 1:length(Names)
  n = Names{i};
  N = str2double(n(2));
  % Stelle den Code-Ordner fest. Siehe auch: serroblib_gen_bitarrays.m
  exp_var = 'S(\d)([RP]+)(\d+)V(\d+)'; % Format "S3RRR1V1" für Varianten
  [tokens_var, ~] = regexp(n,exp_var,'tokens','match');
  if ~isempty(tokens_var) % Name entspricht Variante
    name_haupt = ['S',tokens_var{1}{1},tokens_var{1}{2},tokens_var{1}{3}];
    fcn_dir_var = fullfile(repopath, sprintf('mdl_%ddof', N), name_haupt, ...
                            sprintf('hd_V%s', tokens_var{1}{4}));
    fcn_dir_gen = fullfile(repopath, sprintf('mdl_%ddof', N), name_haupt, ...
                            'hd');
    if exist(fcn_dir_var, 'file')
      fcn_dir = fcn_dir_var;
    else
      fcn_dir = fcn_dir_gen;
    end
  else % Name entspricht Hauptmodell
    fcn_dir = fullfile(repopath, sprintf('mdl_%ddof', N), n, 'hd');
  end
  if ~exist(fcn_dir, 'file')
    warning('%s: Code soll kompiliert werden, obwohl Verzeichnis nicht existiert (%s)', n, fcn_dir);
    continue
  end
  cd(fcn_dir); % wechsele in Verzeichnis (damit Funktionen fürs kompilieren gefunden werden)
  mex_all_matlabfcn_in_dir(fcn_dir, 0);
end

