% Calculate homogenous joint transformation matrices for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6RPRPRR9_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:07:59
% EndTime: 2018-11-23 16:07:59
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->12), mult. (12->12), div. (0->0), fcn. (44->16), ass. (0->17)
t128 = cos(qJ(1));
t127 = cos(qJ(3));
t126 = cos(qJ(5));
t125 = cos(qJ(6));
t124 = sin(qJ(1));
t123 = sin(qJ(3));
t122 = sin(qJ(5));
t121 = sin(qJ(6));
t120 = cos(pkin(6));
t119 = cos(pkin(7));
t118 = cos(pkin(12));
t117 = cos(pkin(13));
t116 = sin(pkin(6));
t115 = sin(pkin(7));
t114 = sin(pkin(12));
t113 = sin(pkin(13));
t1 = [t128, -t124, 0, 0; t124, t128, 0, 0; 0, 0, 1, pkin(8); 0, 0, 0, 1; t118, -t114, 0, pkin(1); t120 * t114, t120 * t118, -t116, -t116 * qJ(2); t116 * t114, t116 * t118, t120, t120 * qJ(2); 0, 0, 0, 1; t127, -t123, 0, pkin(2); t119 * t123, t119 * t127, -t115, -t115 * pkin(9); t115 * t123, t115 * t127, t119, t119 * pkin(9); 0, 0, 0, 1; t117, -t113, 0, pkin(3); t113, t117, 0, 0; 0, 0, 1, qJ(4); 0, 0, 0, 1; t126, -t122, 0, pkin(4); 0, 0, -1, -pkin(10); t122, t126, 0, 0; 0, 0, 0, 1; t125, -t121, 0, pkin(5); 0, 0, -1, -pkin(11); t121, t125, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
