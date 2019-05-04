% Generiere Matlab-Code mit Maple-Dynamik-Toolbox für eine Roboterstruktur
% 
% Eingabe:
% Names
%   Cell-Array mit Liste der Namen der Roboterstrukturen, für die der Code
%   erzeugt werden soll. Der Name entspricht dem Schema "SxRRRyyy" mit
%   x=Anzahl Gelenk-FG und yyy laufende Nummer für "RRR".
% force [1x1 logical]
%   Erzwinge die Neu-Generierung des Codes
% nocopy [1x1 logical]
%   Nur Code generieren, aber nicht zurück kopieren
% mode
%   Modus, welche Dateien generiert werden sollen, damit nicht die
%   vollständige Dynamik jedes mal neu generiert werden muss
%   1: Vollständige Generierung, alles
%   2: Nur aus Vorlagen generierte Funktionen (z.B. Jacobi, inverse Kinematik)
%   3: Alles generieren in Maple, aber keinen Matlab-Code exportieren
%      (wichtig für PKM. Dort wird die Bein-Dynamik nur als Maple gebraucht)
% 
% Vorher: 
% * Funktion maplerepo_path.m muss vorliegen mit Rückgabe des
%   Repo-Pfades der Maple-Dynamik-Toolbox ("HybrDyn")
% * Maple-Eingabedaten müssen für die Roboterstruktur mit
%   serroblib_generate_mapleinput.m erzeugt werden

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function serroblib_generate_code(Names, force, nocopy, mode)

if nargin < 2
  force = false;
end
if nargin < 3
  nocopy = false;
end
if nargin < 4
  mode = 1;
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
  mapleinputfile_tb = fullfile(mrp, 'robot_codegen_definitions', 'robot_env');
  copyfile( mapleinputfile, mapleinputfile_tb );
  if mode == 3
    fid = fopen(mapleinputfile_tb, 'a');
    fprintf(fid, 'codegen_act := false:\n');
    fclose(fid);
  end
  
  % Datei mit Shell-Variablen auch kopieren (wird bei vollständigem
  % Durchlauf der Toolbox überschrieben, aber für Einzelaufruf von Skripten
  % benötigt
  if exist([mapleinputfile,'.sh'], 'file')
    copyfile( [mapleinputfile,'.sh'], [mapleinputfile_tb, '.sh'] );
  elseif mode == 2
    error('Generierung der Vorlagen-Funktionen nicht möglich. %s existiert nicht', [mapleinputfile,'.sh']);
  end
  % Code-Erstellung starten
  if mode == 1 || mode == 3
    fprintf('Starte Code-Generierung %d/%d für %s\n', i, length(Names), n);
    system( sprintf('cd %s && ./robot_codegen_start.sh --fixb_only --notest --parallel', mrp) ); %  > /dev/null
  elseif mode == 2
    fprintf('Generiere Matlab-Funktionen aus Vorlagen (%d/%d) für %s\n', i, length(Names), n);
    system( sprintf('cd %s/robot_codegen_scripts && ./create_git_versioninfo.sh', mrp) );
    system( sprintf('cd %s/robot_codegen_scripts && ./robot_codegen_tmpvar_matlab.sh', mrp) );
    system( sprintf('cd %s/robot_codegen_scripts && ./robot_codegen_matlab_num_varpar.sh', mrp) );
  else
    error('Modus nicht definiert');
  end
  % generierten Code zurückkopieren (alle .m-Dateien)
  % TODO: Hier sind noch viele automatisch generierte Dateien dabei, die
  % eigentlich nicht relevant sind (z.B. eulxyz-Floatbase)
  if ~nocopy
    for f = dir(fullfile(outputdir_tb, '*.m'))'
      copyfile(fullfile(outputdir_tb, f.name), fullfile(outputdir_local, f.name));
    end
  end
  % Definitionen des Roboters zurückkopieren. Damit lassen sich später
  % leichter Roboterspezifische Funktionen neu generieren, ohne die Toolbox
  % neu durchlaufen zu lassen
  copyfile( fullfile(mrp, 'robot_codegen_definitions', 'robot_env.sh'), [mapleinputfile, '.sh']);
end