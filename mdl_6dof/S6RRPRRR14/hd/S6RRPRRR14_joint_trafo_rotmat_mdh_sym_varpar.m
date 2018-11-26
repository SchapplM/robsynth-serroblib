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
% Datum: 2018-11-23 17:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
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
% StartTime: 2018-11-23 17:30:29
% EndTime: 2018-11-23 17:30:29
% DurationCPUTime: 0.04s
% Computational Cost: add. (14->14), mult. (18->18), div. (0->0), fcn. (54->18), ass. (0->19)
t145 = cos(qJ(1));
t144 = cos(qJ(2));
t143 = cos(qJ(4));
t142 = cos(qJ(5));
t141 = cos(qJ(6));
t140 = sin(qJ(1));
t139 = sin(qJ(2));
t138 = sin(qJ(4));
t137 = sin(qJ(5));
t136 = sin(qJ(6));
t135 = cos(pkin(6));
t134 = cos(pkin(7));
t133 = cos(pkin(8));
t132 = cos(pkin(14));
t131 = sin(pkin(6));
t130 = sin(pkin(7));
t129 = sin(pkin(8));
t128 = sin(pkin(14));
t1 = [t145, -t140, 0, 0; t140, t145, 0, 0; 0, 0, 1, pkin(9); 0, 0, 0, 1; t144, -t139, 0, pkin(1); t135 * t139, t135 * t144, -t131, -t131 * pkin(10); t131 * t139, t131 * t144, t135, t135 * pkin(10); 0, 0, 0, 1; t132, -t128, 0, pkin(2); t134 * t128, t134 * t132, -t130, -t130 * qJ(3); t130 * t128, t130 * t132, t134, t134 * qJ(3); 0, 0, 0, 1; t143, -t138, 0, pkin(3); t133 * t138, t133 * t143, -t129, -t129 * pkin(11); t129 * t138, t129 * t143, t133, t133 * pkin(11); 0, 0, 0, 1; t142, -t137, 0, pkin(4); 0, 0, -1, -pkin(12); t137, t142, 0, 0; 0, 0, 0, 1; t141, -t136, 0, pkin(5); 0, 0, -1, -pkin(13); t136, t141, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
