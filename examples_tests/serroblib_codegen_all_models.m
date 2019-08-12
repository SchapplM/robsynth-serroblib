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
serroblib_create_robot_csv_all

%% Code für alle Modelle generieren
only_for_existing_code = true;
for N = 1:7
  % Liste zusammenstellen
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo');
  if ~only_for_existing_code
    II = 1:length(l.Names_Ndof);
  else
    I = (l.AdditionalInfo(:,4) == 1); % Code liegt bereits vor
    II = find(I);
  end
  % alle Modelle generieren
  for iFK = II'
    Name = l.Names_Ndof{iFK};
    serroblib_generate_mapleinput({Name})
    serroblib_generate_code({Name}, true)
    % serroblib_addtopath({Name})
  end
end

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
