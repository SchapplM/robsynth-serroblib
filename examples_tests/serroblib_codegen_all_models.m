% Generieren des Matlab-Codes für alle seriellen Robotermodelle
% 
% Nutzen dieses Skripts:
% * Neu-Generierung aller Matlab-Funktionen, z.B. wenn Änderungen in der
%   HybrDyn-Toolbox erfolgt sind.
% * Generieren von Matlab-Funktionen für Roboter, die neu zur
%   Modelldatenbank hinzugefügt wurden

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover
clear
clc
repopath=fileparts(which('serroblib_path_init.m'));
%% Alle Robotermodell initialisieren
serroblib_gen_bitarrays
% serroblib_create_robot_csv_all

%% Code für alle Modelle generieren
generate_for_existing_code = true;
only_general_models = true;
generate_for_missing_code = true;
only_models_with_offset = true;

for N = 2:7
  % Liste zusammenstellen
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo', 'BitArrays_Ndof');
  I_missing_code = (l.AdditionalInfo(:,4) == 0);
  I_existing_code = (l.AdditionalInfo(:,4) == 1);
  I_general_mdl = (l.AdditionalInfo(:,2) == 0);
  % Prüfe, ob für die Gelenkkoordinaten ein Offset vorliegt. Das ist der
  % Fall, wenn die erste Drehachse nicht die z-Achse ist. (aus Umrechnung
  % von DH nach MDH).
  % die letzten drei Bit geben den Offset an (siehe serroblib_bits2csvline)
  Mask_Offset = uint16(bitshift( uint16(bin2dec('111')), 13)); % `dec2bin(Mask_Offset,16)`
  I_offset = false(length(l.Names_Ndof),1);
  for j = 1:length(l.Names_Ndof) % gehe alle Modelle durch
    for i = 1:N % Prüfe für jedes Gelenk
      Res_Offset = bitand(l.BitArrays_Ndof(j,i), Mask_Offset); % `dec2bin(Res_Offset,16)`
      if Res_Offset ~= 0 % Der Index-Wert Null entspricht dem Eintrag "Kein Offset"
        I_offset(j) = true; 
      end
    end
  end
  % Wende Filterkriterien an
  I = false(length(l.Names_Ndof),1);
  if generate_for_existing_code
    I = I | I_existing_code;
  end
  if generate_for_missing_code
    I = I | I_missing_code;
  end
  if only_general_models
    I = I & I_general_mdl;
  end
  if only_models_with_offset
    I = I & I_offset;
  end
  II = find(I);
  % alle Modelle generieren
  fprintf('Generiere Code für %d/%d Modelle mit %d FG\n', length(II), length(I), N);
  for iFK = II'
    Name = l.Names_Ndof{iFK};
    fprintf('Generiere Matlab-Code für %s\n', Name);
    serroblib_generate_mapleinput({Name})
    serroblib_generate_code({Name}, true)
    % serroblib_addtopath({Name})
  end
end

return % Ab jetzt folgen nur noch Beispiele
%% Nur Vorlagen-Funktionen für alle Modelle neu generieren
% Nur für Modelle, für die bereits Code vorliegt. Sonst wird für alle
% Varianten ohne Code-Ordner auch generiert
for N = 1:7
  % Liste zusammenstellen
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo');
  I = (l.AdditionalInfo(:,4) == 1); % Code liegt bereits vor
  II = find(I);
  % alle Modelle generieren
  j = 0;
  for iFK = II'
    j = j + 1;
    Name = l.Names_Ndof{iFK};
    fprintf('%d/%d (Nr. %d): Generiere Vorlagen für %s\n', j, length(II), iFK, Name);
    % serroblib_generate_code({Name}, true, false, 2)
    serroblib_create_template_functions({Name}, false, true)
  end
end
% Einfacherer Aufruf:
% serroblib_create_template_functions({}, false, true)

%% Code für ein bestimmtes Modell neu generieren
% Beispiel: S6RRPRRR14 (z.B. UPS-Kette)
Name = 'S6RRPRRR14';
serroblib_generate_mapleinput({Name})
serroblib_generate_code({Name}, true)

%% Vorlagen-Funktionen für ein bestimmtes Modell neu generieren
Name = 'S4RRPR1';
% serroblib_generate_code({Name}, true, false, 2)
serroblib_create_template_functions({Name}, false, true)
%% Testfunktionen für ein bestimmtes Modell generieren und ausführen
Name = 'S6RRRRRR10V3';
serroblib_generate_code({Name}, true, false, 4)
