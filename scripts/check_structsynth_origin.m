% Prüfe die Herkunft aller Modelle in der Datenbank.
% Da für jede Kette eine Herkunft angegeben werden kann, darf diese
% Eigenschaft nicht ungesetzt sein.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Benutzereingaben
only_general_model = false; % Prüfe alle Modelle, auch Varianten
%% Initialisierung
roblibpath=fileparts(which('serroblib_path_init.m'));
serroblib_gen_bitarrays(1:7);

%% Alle Modelle in Datenbank durchgehen
for N = 1:7
  fprintf('Prüfe Strukturen mit %d FG\n', N);
  % Alle Roboter aus Datenbank laden
  mdllistfile_Ndof = fullfile(roblibpath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo', 'BitArrays_Origin');
  I_novar = (l.AdditionalInfo(:,2) == 0);
  if only_general_model
    I = find(I_novar);
  else
    I = 1:length(I_novar);
  end
  for j = I(:)'
    Name = l.Names_Ndof{j};
    if l.BitArrays_Origin(j) == 0
      fprintf('Modell %s hat keine Herkunftsangabe\n', Name);
    end
  end
end