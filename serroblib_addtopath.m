% Füge die Modelle gegebener Roboter zum Pfad hinzu

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function serroblib_addtopath(Names)
repopath=fileparts(which('serroblib_path_init.m'));
for i = 1:length(Names)
  n = Names{i};
  N = str2double(n(2)); % Nehme Namensschema SxRRRR mit x=Anzahl FG

  fcn_dir = fullfile(repopath, sprintf('mdl_%ddof', N), n, 'hd');
  addpath(fcn_dir);
end