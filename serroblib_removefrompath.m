% Entferne die Modelle gegebener Roboter vom Matlab-Pfad

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-10
% (C) Institut für mechatronische Systeme, Universität Hannover

function serroblib_removefrompath(Names)
repopath=fileparts(which('serroblib_path_init.m'));
for i = 1:length(Names)
  n = Names{i};
  N = str2double(n(2)); % Nehme Namensschema SxRRRR mit x=Anzahl FG

  fcn_dir = fullfile(repopath, sprintf('mdl_%ddof', N), n, 'hd');
  rmpath(fcn_dir);
end