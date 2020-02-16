% CAD-Modell für Roboter Kuka LBR 4+ vorbereiten
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
% Quelle des CAD-Modells: https://github.com/CentroEPiaggio/kuka-lwr.git

% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-11
% (c) Institut für Regelungstechnik, Universität Hannover

function RS = S7RRRRRRR1_LWR4P_init_CAD(RS, Name, RobName)
N = RS.NJ;
% Pfad zu diesem Ordner bestimmen (wo die CAD-Dateien drin sind)
repopath = fileparts(which('serroblib_path_init.m'));
CAD_basepath = fullfile(repopath, sprintf('mdl_%ddof', N), Name, sprintf('CAD_%s',RobName));
CAD_filenames = {'base.STL', 'link_1.STL', 'link_2.STL', 'link_3.STL', 'link_4.STL', 'link_5.STL', 'link_6.STL', 'link_7.STL'};
CAD_urlpath = 'https://github.com/CentroEPiaggio/kuka-lwr/raw/master/lwr_description/meshes/lwr4plus/visual/';

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
T_mdh_urdf = NaN(4,4,8);
T_mdh_urdf(:,:,1) = eye(4);
T_mdh_urdf(:,:,2) = transl([0;0;-0.200]);
T_mdh_urdf(:,:,3) = trotx(-pi/2);
T_mdh_urdf(:,:,4) = transl([0;0;-0.200]);
T_mdh_urdf(:,:,5) = trotx(pi/2);
T_mdh_urdf(:,:,6) = transl([0;0;-0.195]);
T_mdh_urdf(:,:,7) = trotx(-pi/2);
T_mdh_urdf(:,:,8) = eye(4);%transl([0;0;0.078]);

% Transformation von URDF nach visual KS
% Siehe kuka-lwr/lwr_description/model/kuka_lwr.urdf.xacro
% (Transformation bei visual)
T_urdf_visual = NaN(4,4,8);
T_urdf_visual(:,:,1) = trotz(pi);
T_urdf_visual(:,:,2) = transl([0;0;-0.008])*trotz(pi);
T_urdf_visual(:,:,3) = trotz(pi); % <link name="${name}_2_link">
T_urdf_visual(:,:,4) = transl([0;0;-0.008])*trotz(pi);
T_urdf_visual(:,:,5) = trotz(pi); % <link name="${name}_4_link">
T_urdf_visual(:,:,6) = transl([0;0;-0.008])*trotz(pi);
T_urdf_visual(:,:,7) = trotz(pi);
T_urdf_visual(:,:,8) = trotz(pi);% <link name="${name}_7_link">

for i = 0:7
  if isempty(CAD_filenames{i+1})
    continue
  end
  T_mdh_CAD = T_mdh_urdf(:,:,i+1)*T_urdf_visual(:,:,i+1);
  RS.CAD_add(fullfile(CAD_basepath, CAD_filenames{i+1}), i, T_mdh_CAD);
end