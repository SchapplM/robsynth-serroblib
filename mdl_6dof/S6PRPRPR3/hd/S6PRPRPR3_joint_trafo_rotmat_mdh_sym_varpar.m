% Calculate homogenous joint transformation matrices for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_mdh = S6PRPRPR3_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:56:19
% EndTime: 2018-11-23 14:56:19
% DurationCPUTime: 0.04s
% Computational Cost: add. (10->10), mult. (6->6), div. (0->0), fcn. (30->12), ass. (0->13)
t95 = cos(qJ(2));
t94 = cos(qJ(4));
t93 = cos(qJ(6));
t92 = sin(qJ(2));
t91 = sin(qJ(4));
t90 = sin(qJ(6));
t89 = cos(pkin(6));
t88 = cos(pkin(10));
t87 = cos(pkin(11));
t86 = sin(pkin(6));
t85 = sin(pkin(10));
t84 = sin(pkin(11));
t1 = [t88, -t85, 0, 0; t85, t88, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t95, -t92, 0, pkin(1); t89 * t92, t89 * t95, -t86, -t86 * pkin(7); t86 * t92, t86 * t95, t89, t89 * pkin(7); 0, 0, 0, 1; t87, -t84, 0, pkin(2); t84, t87, 0, 0; 0, 0, 1, qJ(3); 0, 0, 0, 1; t94, -t91, 0, pkin(3); 0, 0, -1, -pkin(8); t91, t94, 0, 0; 0, 0, 0, 1; 0, -1, 0, pkin(4); 0, 0, -1, -qJ(5); 1, 0, 0, 0; 0, 0, 0, 1; t93, -t90, 0, pkin(5); 0, 0, -1, -pkin(9); t90, t93, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
