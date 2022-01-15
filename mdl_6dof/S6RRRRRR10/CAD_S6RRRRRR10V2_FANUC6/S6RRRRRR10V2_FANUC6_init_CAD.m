% CAD-Modell für FANUC6 vorbereiten
% Die Daten liegen im CAD-Repo
%
% Eingabe:
% RS
%   Instanz der Roboterklasse, initialisiert für Roboter
% 
% Ausgabe:
% RS
%   Roboterklasse mit initialisierten Pfaden zum CAD-Modell

% Tobias Blum, MA bei moritz.schappler@imes.uni-hannover.de, 2020-09
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function RS = S6RRRRRR10V2_FANUC6_init_CAD(RS, ~, ~)
N = RS.NJ;
% Pfad zu dem Ordner bestimmen, in dem die CAD-Dateien drin sind
repopath = fileparts(which('robmdlbib_cad_path_init.m'));
if isempty(repopath)
  warning('Das Repo robmdlbib_cad ist nicht im Matlab-Pfad. CAD-Modell nicht verfügbar');
end
CAD_basepath = fullfile(repopath, 'FANUC', 'S6RRRRRR10V2_FANUC6');
if ~exist(CAD_basepath, 'file')
  warning('CAD-Modell liegt nicht im CAD-Repo. Falsche Version?');
end

CAD_filenames = {'base_new.stl', 'segment1_new.stl', 'segment2_new.stl', ...
  'segment3_new.stl', 'segment4_new.stl', 'segment5_new.stl', 'segment6_new.stl'};

% Das Robotermodell wurde mit Körper-Koordinatensystemen nach URDF-Notation
% aufgestellt, die Transformationen der STL-Dateien beziehen sich darauf
% Transformation von MDH-KS nach URDF-KS. Händisch ausprobiert.
T_mdh_urdf = NaN(4,4,7);
T_mdh_urdf(:,:,1) = transl([0;0;-0.5]); % TODO: Nur grob abgeschätzt.
T_mdh_urdf(:,:,2) = eye(4);
T_mdh_urdf(:,:,3) = transl([0.870;0;-0.121-0.01]);
T_mdh_urdf(:,:,4) = trotx(pi)*trotz(pi);
T_mdh_urdf(:,:,5) = trotz(pi/2); 
T_mdh_urdf(:,:,6) = trotz(pi); 
T_mdh_urdf(:,:,7) = eye(4);  

for i = 0:6
  if isempty(CAD_filenames{i+1})
    continue
  end
  T_mdh_CAD = T_mdh_urdf(:,:,i+1);
  RS.CAD_add(fullfile(CAD_basepath, CAD_filenames{i+1}), i, T_mdh_CAD);
end