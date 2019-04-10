% Generiere die csv-Dateien für alle Robotermodelle (zum Abspeichern der
% Strukturen mit konkreten Zahlenwerten für die Parameter)
% 
% Schreibt Dateien:
% /mdl_xdof/SRR...PRy/models.csv (für alle x FG, RRR... Gelenkanordnungen
% und y Modellvariationen

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für Mechatronische Systeme, Universität Hannover

repopath=fileparts(which('serroblib_path_init.m'));

for N = 1:7
  % Bit-Arrays aktualisieren
  serroblib_gen_bitarrays(N);
  
  fprintf('Erstelle Modelldateien für Roboter mit %d FG\n', N);
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'AdditionalInfo');
  for j = 1:length(l.Names_Ndof)
    isvariant = l.AdditionalInfo(j,2);
    variantof = l.AdditionalInfo(j,3);
    if isvariant
      fprintf('Modell %s ist nur eine Variante von %s. Nichts erstellen.\n', ...
        l.Names_Ndof{j}, l.Names_Ndof{variantof});
      continue
    end
    [exists, filepath_csv] = serroblib_create_robot_csv(l.Names_Ndof{j});
    if exists
      fprintf('Datei %s existiert bereits\n', filepath_csv);
    else
      fprintf('Datei %s erstellt\n', filepath_csv);
    end
  end
end