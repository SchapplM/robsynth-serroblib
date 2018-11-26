% Calculate homogenous joint transformation matrices for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6RPRRRR12_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:40:40
% EndTime: 2018-11-23 16:40:40
% DurationCPUTime: 0.05s
% Computational Cost: add. (14->14), mult. (18->18), div. (0->0), fcn. (54->18), ass. (0->19)
t147 = cos(qJ(1));
t146 = cos(qJ(3));
t145 = cos(qJ(4));
t144 = cos(qJ(5));
t143 = cos(qJ(6));
t142 = sin(qJ(1));
t141 = sin(qJ(3));
t140 = sin(qJ(4));
t139 = sin(qJ(5));
t138 = sin(qJ(6));
t137 = cos(pkin(6));
t136 = cos(pkin(7));
t135 = cos(pkin(8));
t134 = cos(pkin(14));
t133 = sin(pkin(6));
t132 = sin(pkin(7));
t131 = sin(pkin(8));
t130 = sin(pkin(14));
t1 = [t147, -t142, 0, 0; t142, t147, 0, 0; 0, 0, 1, pkin(9); 0, 0, 0, 1; t134, -t130, 0, pkin(1); t137 * t130, t137 * t134, -t133, -t133 * qJ(2); t133 * t130, t133 * t134, t137, t137 * qJ(2); 0, 0, 0, 1; t146, -t141, 0, pkin(2); t136 * t141, t136 * t146, -t132, -t132 * pkin(10); t132 * t141, t132 * t146, t136, t136 * pkin(10); 0, 0, 0, 1; t145, -t140, 0, pkin(3); t135 * t140, t135 * t145, -t131, -t131 * pkin(11); t131 * t140, t131 * t145, t135, t135 * pkin(11); 0, 0, 0, 1; t144, -t139, 0, pkin(4); 0, 0, -1, -pkin(12); t139, t144, 0, 0; 0, 0, 0, 1; t143, -t138, 0, pkin(5); 0, 0, -1, -pkin(13); t138, t143, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
