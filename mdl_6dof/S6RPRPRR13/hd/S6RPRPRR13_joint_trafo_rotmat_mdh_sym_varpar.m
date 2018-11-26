% Calculate homogenous joint transformation matrices for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6RPRPRR13_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:10:25
% EndTime: 2018-11-23 16:10:25
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->12), mult. (12->12), div. (0->0), fcn. (40->14), ass. (0->15)
t113 = cos(qJ(1));
t112 = cos(qJ(3));
t111 = cos(qJ(5));
t110 = cos(qJ(6));
t109 = sin(qJ(1));
t108 = sin(qJ(3));
t107 = sin(qJ(5));
t106 = sin(qJ(6));
t105 = cos(pkin(6));
t104 = cos(pkin(7));
t103 = cos(pkin(12));
t102 = sin(pkin(6));
t101 = sin(pkin(7));
t100 = sin(pkin(12));
t1 = [t113, -t109, 0, 0; t109, t113, 0, 0; 0, 0, 1, pkin(8); 0, 0, 0, 1; t103, -t100, 0, pkin(1); t105 * t100, t105 * t103, -t102, -t102 * qJ(2); t102 * t100, t102 * t103, t105, t105 * qJ(2); 0, 0, 0, 1; t112, -t108, 0, pkin(2); t104 * t108, t104 * t112, -t101, -t101 * pkin(9); t101 * t108, t101 * t112, t104, t104 * pkin(9); 0, 0, 0, 1; 0, -1, 0, pkin(3); 0, 0, -1, -qJ(4); 1, 0, 0, 0; 0, 0, 0, 1; t111, -t107, 0, pkin(4); 0, 0, -1, -pkin(10); t107, t111, 0, 0; 0, 0, 0, 1; t110, -t106, 0, pkin(5); 0, 0, -1, -pkin(11); t106, t110, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
