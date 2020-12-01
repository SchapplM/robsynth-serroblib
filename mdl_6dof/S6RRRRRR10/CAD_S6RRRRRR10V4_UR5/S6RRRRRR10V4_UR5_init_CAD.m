% CAD-Modell für UR5 vorbereiten
%
% Eingabe:
% RS
%   Instanz der Roboterklasse, initialisiert für Roboter
% Name
%   Name des Robotermodells nach dem Schema "SxRRRyy"
% RobName
%   Name der Roboterparameter. Das CAD-Modell bezieht sich auf die
%   Roboterparameter mit dem gleichen Namen
% 
% Ausgabe:
% RS
%   Roboterklasse mit initialisierten Pfaden zum CAD-Modell
%
% Siehe auch: serroblib_create_robot_class.m
% Quelle des CAD-Modells: https://github.com/ros-industrial/universal_robot

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-09
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function RS = S6RRRRRR10V4_UR5_init_CAD(RS, Name, RobName)
N = RS.NJ;
% Pfad zu diesem Ordner bestimmen (wo die CAD-Dateien drin sind)
repopath = fileparts(which('serroblib_path_init.m'));
CAD_basepath = fullfile(repopath, sprintf('mdl_%ddof', N), Name, sprintf('CAD_%s',RobName));
CAD_filenames = {'base.stl', 'shoulder.stl', 'upperarm.stl', 'forearm.stl', 'wrist1.stl', 'wrist2.stl', 'wrist3.stl'};
CAD_urlpath = 'https://github.com/ros-industrial/universal_robot/raw/kinetic-devel/ur_description/meshes/ur5/collision/';

% CAD-Dateien herunterladen. Die Dateien sollen nicht in diesem Repo
% eingecheckt werden, damit es nicht zu Lizenproblemen kommt
for i = 1:length(CAD_filenames)
  if ~exist(fullfile(CAD_basepath, CAD_filenames{i}), 'file')
    try  % Versuche, herunterzuladen, falls Datei fehlt
      websave(fullfile(CAD_basepath, CAD_filenames{i}), ...
              [CAD_urlpath,  CAD_filenames{i}]); % URL to file on GitHub
      fprintf('Datei %s von Github heruntergeladen\n', CAD_filenames{i});
    catch
      warning('Fehler beim Herunterladen der Datei: %s', CAD_filenames{i});
      CAD_filenames{i} = '';
    end
  end
end

% Das Robotermodell wurde mit Körper-Koordinatensystemen nach URDF-Notation
% aufgestellt, die Transformationen der STL-Dateien beziehen sich darauf
% Transformation von MDH-KS nach URDF-KS. Händisch ausprobiert.
T_mdh_urdf = NaN(4,4,7);
T_mdh_urdf(:,:,1) = eye(4);
T_mdh_urdf(:,:,2) = trotx(-pi/2);
T_mdh_urdf(:,:,3) = trotz(-pi/2)*transl([0;0;0.1304]);
T_mdh_urdf(:,:,4) = trotx(pi/2)*transl([0;0.0107;0]); 
T_mdh_urdf(:,:,5) = trotx(pi/2)*transl(-[0;0.0984;0]); 
T_mdh_urdf(:,:,6) = eye(4)*transl([0;0;-0.0963]); 
T_mdh_urdf(:,:,7) = trotx(-pi/2)*transl([0;-0.07697;0]);

% Transformation von URDF nach visual KS
% Siehe ur_description/urdf/ur5.urdf.xacro auf GitHub (Transformation bei visual)
T_urdf_visual = NaN(4,4,7);
T_urdf_visual(:,:,1) = eye(4);
T_urdf_visual(:,:,2) = eye(4);
T_urdf_visual(:,:,3) = trotx(pi/2);
T_urdf_visual(:,:,4) = troty(-pi/2);
T_urdf_visual(:,:,5) = eye(4);
T_urdf_visual(:,:,6) = eye(4);
T_urdf_visual(:,:,7) = eye(4);
T_urdf_visual(:,:,8) = eye(4);
% In Roboterklasse eintragen
for i = 0:6
  if isempty(CAD_filenames{i+1})
    continue
  end
  T_mdh_CAD = T_mdh_urdf(:,:,i+1)*T_urdf_visual(:,:,i+1);
  RS.CAD_add(fullfile(CAD_basepath, CAD_filenames{i+1}), i, T_mdh_CAD);
end