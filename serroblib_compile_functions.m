% Kompiliere alle Matlab-Funktionen aller Robotermodelle in der Bibliothek
% 
% Eingabe:
% Names
%   Cell array mit Namen der zu kompilierenden Robotermodelle
%   Optional: Wenn nicht angegeben, werden alle kompiliert.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function serroblib_compile_functions(Names)

% Prüfe Eingabeargument
repopath=fileparts(which('serroblib_path_init.m'));
if nargin == 0
  % Stelle Liste aller Roboter zusammen
  Names = {};
  for N = 1:6
    mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
    if ~exist(mdllistfile_Ndof, 'file')
      continue
    end
    % Lade in .mat gespeicherte Datenbank
    l = load(mdllistfile_Ndof, 'Names_Ndof');
    Names = {Names{:}, l.Names_Ndof{:}};
  end
end

% Gehe alle Modellnamen durch und kompiliere alle m-Dateien in dem
% jeweiligen Ordner.
for i = 1:length(Names)
  n = Names{i};
  N = str2double(n(2));
 
  fcn_dir = fullfile(repopath, sprintf('mdl_%ddof', N), n, 'hd');
  cd(fcn_dir); % wechsele in Verzeichnis (damit Funktionen fürs kompilieren gefunden werden)
  mex_all_matlabfcn_in_dir(fcn_dir, 0);
end

