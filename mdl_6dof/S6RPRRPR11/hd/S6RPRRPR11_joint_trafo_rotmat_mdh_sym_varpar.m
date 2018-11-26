% Calculate homogenous joint transformation matrices for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6RPRRPR11_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:22:04
% EndTime: 2018-11-23 16:22:04
% DurationCPUTime: 0.05s
% Computational Cost: add. (12->12), mult. (12->12), div. (0->0), fcn. (44->16), ass. (0->17)
t116 = cos(qJ(1));
t115 = cos(qJ(3));
t114 = cos(qJ(4));
t113 = cos(qJ(6));
t112 = sin(qJ(1));
t111 = sin(qJ(3));
t110 = sin(qJ(4));
t109 = sin(qJ(6));
t108 = cos(pkin(6));
t107 = cos(pkin(7));
t106 = cos(pkin(12));
t105 = cos(pkin(13));
t104 = sin(pkin(6));
t103 = sin(pkin(7));
t102 = sin(pkin(12));
t101 = sin(pkin(13));
t1 = [t116, -t112, 0, 0; t112, t116, 0, 0; 0, 0, 1, pkin(8); 0, 0, 0, 1; t106, -t102, 0, pkin(1); t108 * t102, t108 * t106, -t104, -t104 * qJ(2); t104 * t102, t104 * t106, t108, t108 * qJ(2); 0, 0, 0, 1; t115, -t111, 0, pkin(2); t107 * t111, t107 * t115, -t103, -t103 * pkin(9); t103 * t111, t103 * t115, t107, t107 * pkin(9); 0, 0, 0, 1; t114, -t110, 0, pkin(3); 0, 0, -1, -pkin(10); t110, t114, 0, 0; 0, 0, 0, 1; t105, -t101, 0, pkin(4); 0, 0, -1, -qJ(5); t101, t105, 0, 0; 0, 0, 0, 1; t113, -t109, 0, pkin(5); t109, t113, 0, 0; 0, 0, 1, pkin(11); 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
