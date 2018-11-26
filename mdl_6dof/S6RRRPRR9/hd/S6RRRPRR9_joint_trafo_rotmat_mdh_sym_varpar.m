% Calculate homogenous joint transformation matrices for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6RRRPRR9_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:56:55
% EndTime: 2018-11-23 17:56:55
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->12), mult. (12->12), div. (0->0), fcn. (44->16), ass. (0->17)
t127 = cos(qJ(1));
t126 = cos(qJ(2));
t125 = cos(qJ(3));
t124 = cos(qJ(5));
t123 = cos(qJ(6));
t122 = sin(qJ(1));
t121 = sin(qJ(2));
t120 = sin(qJ(3));
t119 = sin(qJ(5));
t118 = sin(qJ(6));
t117 = cos(pkin(6));
t116 = cos(pkin(7));
t115 = cos(pkin(13));
t114 = sin(pkin(6));
t113 = sin(pkin(7));
t112 = sin(pkin(13));
t1 = [t127, -t122, 0, 0; t122, t127, 0, 0; 0, 0, 1, pkin(8); 0, 0, 0, 1; t126, -t121, 0, pkin(1); t117 * t121, t117 * t126, -t114, -t114 * pkin(9); t114 * t121, t114 * t126, t117, t117 * pkin(9); 0, 0, 0, 1; t125, -t120, 0, pkin(2); t116 * t120, t116 * t125, -t113, -t113 * pkin(10); t113 * t120, t113 * t125, t116, t116 * pkin(10); 0, 0, 0, 1; t115, -t112, 0, pkin(3); t112, t115, 0, 0; 0, 0, 1, qJ(4); 0, 0, 0, 1; t124, -t119, 0, pkin(4); 0, 0, -1, -pkin(11); t119, t124, 0, 0; 0, 0, 0, 1; t123, -t118, 0, pkin(5); 0, 0, -1, -pkin(12); t118, t123, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
