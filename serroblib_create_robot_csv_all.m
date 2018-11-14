% Generiere die csv-Dateien für alle Robotermodelle (zum Abspeichern der
% Strukturen mit konkreten Zahlenwerten für die Parameter)
% 
% Schreibt Dateien:
% /mdl_xdof/SRR...PRy/models.csv (für alle x FG, RRR... Gelenkanordnungen
% und y Modellvariationen

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

for N = 1:6
  % Bit-Arrays aktualisieren
  serroblib_gen_bitarrays(N);
  
  fprintf('Erstelle Modelldateien für Roboter mit %d FG\n', N);
  repopath=fileparts(which('serroblib_path_init.m'));
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof');
  for j = 1:length(l.Names_Ndof)
    [exists, filepath_csv] = serroblib_create_robot_csv(l.Names_Ndof{j});
    if exists
      fprintf('Datei %s existiert bereits\n', filepath_csv);
    else
      fprintf('Datei %s erstellt\n', filepath_csv);
    end
  end
end