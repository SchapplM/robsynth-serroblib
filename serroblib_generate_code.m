% Generiere Matlab-Code mit Maple-Dynamik-Toolbox für eine Roboterstruktur
% 
% Eingabe:
% Names
%   Cell-Array mit Liste der Namen der Roboterstrukturen, für die der Code
%   erzeugt werden soll. Der Name entspricht dem Schema "SxRRRyyy" mit
%   x=Anzahl Gelenk-FG und yyy laufende Nummer für "RRR".
% force [1x1 logical]
%   Erzwinge die Neu-Generierung des Codes
% 
% Vorher: 
% * Funktion maplerepo_path.m muss vorliegen mit Rückgabe des
%   Repo-Pfades der Maple-Dynamik-Toolbox ("HybrDyn")
% * Maple-Eingabedaten müssen für die Roboterstruktur mit
%   serroblib_generate_mapleinput.m erzeugt werden

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function serroblib_generate_code(Names, force)

if nargin < 2
  force = false;
end

repopath=fileparts(which('serroblib_path_init.m'));

for i = 1:length(Names)
  n = Names{i};
  N = str2double(n(2));

  % Pfad zur Maple-Dynamik-Toolbox (muss im Repo abgelegt werden)
  mrp = maplerepo_path();
  
  % Maple-Toolbox-Eingabe laden (wurde an anderer Stelle erzeugt)
  % (durch serroblib_generate_mapleinput.m)
  mapleinputfile=fullfile(repopath, sprintf('mdl_%ddof', N), n, 'hd', sprintf('robot_env_%s', n));
  if ~exist(mapleinputfile, 'file')
    error('Datei %s existiert nicht. Wurde `serroblib_generate_mapleinput.m` ausgeführt?', fileparts(mapleinputfile) );
  end
  % Verzeichnisse für die zu erzeugenden Matlab-Funktionen
  outputdir_tb = fullfile(mrp, 'codeexport', n, 'matlabfcn'); % Verzeichnis in der Maple-Toolbox
  outputdir_local = fileparts(mapleinputfile); % Verzeichnis in der Bibliothek
  
  % Prüfe, ob Code schon einmal generiert wurde 
  % (und im Zielverzeichnis vorliegt)
  if ~force && length(dir(fullfile(outputdir_local, '*.m'))) > 20
    % das werden wohl schon genug .m-Dateien sein.
    continue
  end
  
  % Eingabedatei kopieren
  copyfile( mapleinputfile, fullfile(mrp, 'robot_codegen_definitions', 'robot_env') );
  
  % Code-Erstellung starten
  fprintf('Starte Code-Generierung %d/%d für %s\n', i, length(Names), n);
  system( sprintf('cd %s && ./robot_codegen_start.sh --fixb_only --notest --parallel', mrp) ); %  > /dev/null
  
  % generierten Code zurückkopieren (alle .m-Dateien)
  for f = dir(fullfile(outputdir_tb, '*.m'))'
    copyfile(fullfile(outputdir_tb, f.name), fullfile(outputdir_local, f.name));
  end
end