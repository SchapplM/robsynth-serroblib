% Calculate homogenous joint transformation matrices for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6PPRRPR1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:49:43
% EndTime: 2018-11-23 14:49:43
% DurationCPUTime: 0.05s
% Computational Cost: add. (12->12), mult. (12->12), div. (0->0), fcn. (44->16), ass. (0->17)
t116 = cos(qJ(3));
t115 = cos(qJ(4));
t114 = cos(qJ(6));
t113 = sin(qJ(3));
t112 = sin(qJ(4));
t111 = sin(qJ(6));
t110 = cos(pkin(6));
t109 = cos(pkin(7));
t108 = cos(pkin(11));
t107 = cos(pkin(12));
t106 = cos(pkin(13));
t105 = sin(pkin(6));
t104 = sin(pkin(7));
t103 = sin(pkin(11));
t102 = sin(pkin(12));
t101 = sin(pkin(13));
t1 = [t108, -t103, 0, 0; t103, t108, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t107, -t102, 0, pkin(1); t110 * t102, t110 * t107, -t105, -t105 * qJ(2); t105 * t102, t105 * t107, t110, t110 * qJ(2); 0, 0, 0, 1; t116, -t113, 0, pkin(2); t109 * t113, t109 * t116, -t104, -t104 * pkin(8); t104 * t113, t104 * t116, t109, t109 * pkin(8); 0, 0, 0, 1; t115, -t112, 0, pkin(3); 0, 0, -1, -pkin(9); t112, t115, 0, 0; 0, 0, 0, 1; t106, -t101, 0, pkin(4); 0, 0, -1, -qJ(5); t101, t106, 0, 0; 0, 0, 0, 1; t114, -t111, 0, pkin(5); t111, t114, 0, 0; 0, 0, 1, pkin(10); 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
