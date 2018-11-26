% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPPRRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:50:59
% EndTime: 2018-11-23 15:50:59
% DurationCPUTime: 0.11s
% Computational Cost: add. (142->50), mult. (87->46), div. (0->0), fcn. (132->10), ass. (0->35)
t16 = sin(pkin(10));
t40 = t16 * pkin(3);
t20 = sin(qJ(1));
t14 = pkin(10) + qJ(4);
t8 = qJ(5) + t14;
t5 = cos(t8);
t39 = t20 * t5;
t22 = cos(qJ(1));
t38 = t22 * t5;
t37 = t20 * t16;
t19 = sin(qJ(6));
t36 = t20 * t19;
t21 = cos(qJ(6));
t35 = t20 * t21;
t34 = t22 * t19;
t33 = t22 * t21;
t18 = -pkin(7) - qJ(3);
t6 = sin(t14);
t2 = pkin(4) * t6 + t40;
t32 = -qJ(2) - t2;
t15 = pkin(6) + 0;
t31 = t20 * pkin(1) + 0;
t30 = pkin(2) + t15;
t29 = t22 * pkin(1) + t20 * qJ(2) + 0;
t17 = cos(pkin(10));
t28 = t17 * pkin(3) + t30;
t4 = sin(t8);
t27 = pkin(5) * t4 - pkin(9) * t5;
t7 = cos(t14);
t26 = pkin(4) * t7 + t28;
t13 = -pkin(8) + t18;
t25 = -t20 * t13 + t31;
t24 = -t22 * qJ(2) + t31;
t23 = -t22 * t13 + t20 * t2 + t29;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t20, 0, 0; t20, t22, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; 0, -t22, t20, t29; 0, -t20, -t22, t24; 1, 0, 0, t15; 0, 0, 0, 1; t37, t20 * t17, t22, t22 * qJ(3) + t29; -t22 * t16, -t22 * t17, t20, t20 * qJ(3) + t24; t17, -t16, 0, t30; 0, 0, 0, 1; t20 * t6, t20 * t7, t22, pkin(3) * t37 - t22 * t18 + t29; -t22 * t6, -t22 * t7, t20, -t20 * t18 + (-qJ(2) - t40) * t22 + t31; t7, -t6, 0, t28; 0, 0, 0, 1; t20 * t4, t39, t22, t23; -t22 * t4, -t38, t20, t32 * t22 + t25; t5, -t4, 0, t26; 0, 0, 0, 1; t4 * t35 + t34, -t4 * t36 + t33, -t39, t27 * t20 + t23; -t4 * t33 + t36, t4 * t34 + t35, t38 (-t27 + t32) * t22 + t25; t5 * t21, -t5 * t19, t4, t5 * pkin(5) + t4 * pkin(9) + t26; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
