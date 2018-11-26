% Calculate homogenous joint transformation matrices for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6PPPRRR1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:48:38
% EndTime: 2018-11-23 14:48:38
% DurationCPUTime: 0.04s
% Computational Cost: add. (14->14), mult. (18->18), div. (0->0), fcn. (54->18), ass. (0->19)
t146 = cos(qJ(4));
t145 = cos(qJ(5));
t144 = cos(qJ(6));
t143 = sin(qJ(4));
t142 = sin(qJ(5));
t141 = sin(qJ(6));
t140 = cos(pkin(6));
t139 = cos(pkin(7));
t138 = cos(pkin(8));
t137 = cos(pkin(12));
t136 = cos(pkin(13));
t135 = cos(pkin(14));
t134 = sin(pkin(6));
t133 = sin(pkin(7));
t132 = sin(pkin(8));
t131 = sin(pkin(12));
t130 = sin(pkin(13));
t129 = sin(pkin(14));
t1 = [t137, -t131, 0, 0; t131, t137, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t136, -t130, 0, pkin(1); t140 * t130, t140 * t136, -t134, -t134 * qJ(2); t134 * t130, t134 * t136, t140, t140 * qJ(2); 0, 0, 0, 1; t135, -t129, 0, pkin(2); t139 * t129, t139 * t135, -t133, -t133 * qJ(3); t133 * t129, t133 * t135, t139, t139 * qJ(3); 0, 0, 0, 1; t146, -t143, 0, pkin(3); t138 * t143, t138 * t146, -t132, -t132 * pkin(9); t132 * t143, t132 * t146, t138, t138 * pkin(9); 0, 0, 0, 1; t145, -t142, 0, pkin(4); 0, 0, -1, -pkin(10); t142, t145, 0, 0; 0, 0, 0, 1; t144, -t141, 0, pkin(5); 0, 0, -1, -pkin(11); t141, t144, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
