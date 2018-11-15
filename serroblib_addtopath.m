% Füge die Modelle gegebener Roboter zum Pfad hinzu
% 
% Eingabe:
% Names
%   Cell-Array mit Namen der Robotermodelle, deren Funktionen zum
%   Matlab-Pfad hinzugefügt werden sollen.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function serroblib_addtopath(Names)
if ~iscell(Names)
  error('serroblib_addtopath: Eingegebene Roboternamen müssen ein Cell-Array sein');
end
repopath=fileparts(which('serroblib_path_init.m'));
for i = 1:length(Names)
  n = Names{i};
  N = str2double(n(2)); % Nehme Namensschema SxRRRR mit x=Anzahl FG

  fcn_dir = fullfile(repopath, sprintf('mdl_%ddof', N), n, 'hd');
  if exist(fcn_dir, 'file')
    addpath(fcn_dir);
  end
end