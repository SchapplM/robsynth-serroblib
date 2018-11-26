% Calculate homogenous joint transformation matrices for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6PPRPRR1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:49:06
% EndTime: 2018-11-23 14:49:06
% DurationCPUTime: 0.05s
% Computational Cost: add. (12->12), mult. (12->12), div. (0->0), fcn. (44->16), ass. (0->17)
t128 = cos(qJ(3));
t127 = cos(qJ(5));
t126 = cos(qJ(6));
t125 = sin(qJ(3));
t124 = sin(qJ(5));
t123 = sin(qJ(6));
t122 = cos(pkin(6));
t121 = cos(pkin(7));
t120 = cos(pkin(11));
t119 = cos(pkin(12));
t118 = cos(pkin(13));
t117 = sin(pkin(6));
t116 = sin(pkin(7));
t115 = sin(pkin(11));
t114 = sin(pkin(12));
t113 = sin(pkin(13));
t1 = [t120, -t115, 0, 0; t115, t120, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t119, -t114, 0, pkin(1); t122 * t114, t122 * t119, -t117, -t117 * qJ(2); t117 * t114, t117 * t119, t122, t122 * qJ(2); 0, 0, 0, 1; t128, -t125, 0, pkin(2); t121 * t125, t121 * t128, -t116, -t116 * pkin(8); t116 * t125, t116 * t128, t121, t121 * pkin(8); 0, 0, 0, 1; t118, -t113, 0, pkin(3); t113, t118, 0, 0; 0, 0, 1, qJ(4); 0, 0, 0, 1; t127, -t124, 0, pkin(4); 0, 0, -1, -pkin(9); t124, t127, 0, 0; 0, 0, 0, 1; t126, -t123, 0, pkin(5); 0, 0, -1, -pkin(10); t123, t126, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
