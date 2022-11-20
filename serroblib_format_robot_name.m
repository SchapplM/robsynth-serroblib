% Formatiere einen Roboternamen für Publikationen (in Latex)
% 
% Eingabe:
% Name:
%   Name des Roboters
% 
% Ausgabe:
% FormatName:
%   Latex-Format, Akzente über Buchstaben für Parallelität der Gelenke
% 
% Siehe: parroblib_format_robot_name.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-11
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function FormatName = serroblib_format_robot_name(RobName)
FormatName = '';

%% Roboter-Klasse initialisieren und Parallelität prüfen
R = serroblib_create_robot_class(RobName);
R.fill_fcn_handles(false, false);
NLegJ = R.NJ;
TSS = R.gen_testsettings(false, true);
% Benutze Funktion aus der Maßsynthese-Toolbox
Set_tmp = struct('general', struct('matfile_verbosity', 0));
Structure_tmp = struct('Number', 1, 'Name', RobName);
pgroups_orig = cds_joint_parallelity(R, Set_tmp, Structure_tmp, TSS.Q);
pgroups = pgroups_orig;
%% Zeichenkette für Parallelität generieren. Siehe Kong/Gosselin 2007, S.10
% Kopie aus parroblib_format_robot_name
Chain_Name = RobName;
% Variablen mit Latex-Code für Roboter-Namen
Chain_StructName = '';
% Setze die Nummer 0 für Gelenke, die zu keiner parallelen Gruppe gehören
for j = 1:NLegJ
  if sum(pgroups == pgroups(j)) == 1
    pgroups(j) = 0; % Dadurch dann kein Akzent auf dem Buchstaben
  end
end
% Entferne nicht belegte Nummern
for k = 1:length(pgroups) % mehrfach, falls größere Lücke
  for j = 1:max(pgroups)
    if ~any(pgroups==j) % reduziere alle folgenden Nummern um 1
      pgroups(pgroups>j) = pgroups(pgroups>j) - 1;
    end
  end
end
% Setze die Markierungen entsprechend der Gruppen
for j = 1:NLegJ
  groupidx = pgroups(j); % hochzählen
  if groupidx == 0
    % Diese Gelenkausrichtung gibt es nur einmal. Es muss kein
    % Gruppensymbol darüber gelegt werden
    newsymbol = '{';
  elseif groupidx == 1
    newsymbol = '{\`';
  elseif groupidx == 2
    newsymbol = '{\''';
  elseif groupidx == 3
    newsymbol = '{\=';
  else
    error('mehr als drei Achsrichtungen nicht vorgesehen');
  end
  % Füge "P"/"R" hinzu
  newsymbol = [newsymbol, Chain_Name(2+j), '}']; %#ok<AGROW>
  Chain_StructName = [Chain_StructName, newsymbol]; %#ok<AGROW>
end
FormatName = Chain_StructName;
%% Ende
if isempty(FormatName)
  error('Modus nicht implementiert');
end
