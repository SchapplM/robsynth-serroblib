% Calculate homogenous joint transformation matrices for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6RRPRRR14_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:07:39
% EndTime: 2019-01-03 10:07:39
% DurationCPUTime: 0.05s
% Computational Cost: add. (14->14), mult. (18->18), div. (0->0), fcn. (54->18), ass. (0->19)
t128 = cos(qJ(1));
t127 = cos(qJ(2));
t126 = cos(qJ(4));
t125 = cos(qJ(5));
t124 = cos(qJ(6));
t123 = sin(qJ(1));
t122 = sin(qJ(2));
t121 = sin(qJ(4));
t120 = sin(qJ(5));
t119 = sin(qJ(6));
t118 = cos(pkin(6));
t117 = cos(pkin(7));
t116 = cos(pkin(8));
t115 = cos(pkin(14));
t114 = sin(pkin(6));
t113 = sin(pkin(7));
t112 = sin(pkin(8));
t111 = sin(pkin(14));
t1 = [t128, -t123, 0, 0; t123, t128, 0, 0; 0, 0, 1, pkin(9); 0, 0, 0, 1; t127, -t122, 0, pkin(1); t118 * t122, t118 * t127, -t114, -t114 * pkin(10); t114 * t122, t114 * t127, t118, t118 * pkin(10); 0, 0, 0, 1; t115, -t111, 0, pkin(2); t117 * t111, t117 * t115, -t113, -t113 * qJ(3); t113 * t111, t113 * t115, t117, t117 * qJ(3); 0, 0, 0, 1; t126, -t121, 0, pkin(3); t116 * t121, t116 * t126, -t112, -t112 * pkin(11); t112 * t121, t112 * t126, t116, t116 * pkin(11); 0, 0, 0, 1; t125, -t120, 0, pkin(4); 0, 0, -1, -pkin(12); t120, t125, 0, 0; 0, 0, 0, 1; t124, -t119, 0, pkin(5); 0, 0, -1, -pkin(13); t119, t124, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
