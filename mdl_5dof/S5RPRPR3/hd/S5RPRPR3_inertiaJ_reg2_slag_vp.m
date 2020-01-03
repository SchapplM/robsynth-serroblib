% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t40 = sin(pkin(9));
t66 = -0.2e1 * t40;
t42 = cos(pkin(9));
t65 = 0.2e1 * t42;
t43 = cos(pkin(8));
t61 = t43 * pkin(1);
t32 = pkin(2) + t61;
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t41 = sin(pkin(8));
t62 = t41 * pkin(1);
t18 = t45 * t32 + t47 * t62;
t15 = qJ(4) + t18;
t17 = t47 * t32 - t45 * t62;
t20 = -t42 * pkin(4) - t40 * pkin(7) - pkin(3);
t4 = -t17 + t20;
t44 = sin(qJ(5));
t46 = cos(qJ(5));
t56 = t46 * t42;
t3 = t15 * t56 + t44 * t4;
t36 = t40 ^ 2;
t58 = t36 * t46;
t64 = t15 * t58 + t3 * t42;
t29 = t44 * t42;
t2 = -t15 * t29 + t46 * t4;
t52 = qJ(4) * t42;
t9 = t46 * t20 - t44 * t52;
t63 = -t2 - t9;
t16 = -pkin(3) - t17;
t60 = pkin(3) - t16;
t10 = t44 * t20 + t46 * t52;
t34 = t36 * qJ(4);
t59 = t10 * t42 + t46 * t34;
t12 = t36 * t15;
t57 = t44 * t40;
t37 = t42 ^ 2;
t55 = t37 * t15 + t12;
t35 = t37 * qJ(4);
t54 = t35 + t34;
t38 = t44 ^ 2;
t39 = t46 ^ 2;
t53 = t38 + t39;
t26 = t40 * t65;
t51 = t2 * t46 + t3 * t44;
t50 = t10 * t44 + t9 * t46;
t48 = qJ(4) ^ 2;
t33 = t36 * t48;
t31 = t46 * t40;
t30 = t39 * t36;
t28 = t38 * t36;
t24 = t44 * t34;
t23 = -0.2e1 * t44 * t58;
t22 = t56 * t66;
t21 = t44 * t26;
t19 = t53 * t40;
t14 = t15 ^ 2;
t11 = t36 * t14;
t8 = t15 * t34;
t6 = t44 * t12;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t61, -0.2e1 * t62, 0, (t41 ^ 2 + t43 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t17, -0.2e1 * t18, 0, t17 ^ 2 + t18 ^ 2, t36, t26, 0, t37, 0, 0, -0.2e1 * t16 * t42, 0.2e1 * t16 * t40, 0.2e1 * t55, t37 * t14 + t16 ^ 2 + t11, t30, t23, t22, t28, t21, t37, -0.2e1 * t2 * t42 + 0.2e1 * t6, 0.2e1 * t64, t51 * t66, t2 ^ 2 + t3 ^ 2 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t15 * t42 - t2 * t44 + t3 * t46) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 + t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 + t28 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t17, -t18, 0, 0, t36, t26, 0, t37, 0, 0, t60 * t42, -t60 * t40, t54 + t55, -t16 * pkin(3) + t15 * t35 + t8, t30, t23, t22, t28, t21, t37, t63 * t42 + t24 + t6, t59 + t64, (t63 * t46 + (-t10 - t3) * t44) * t40, t3 * t10 + t2 * t9 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t10 * t46 - t44 * t9 - t52) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, t26, 0, t37, 0, 0, pkin(3) * t65, pkin(3) * t66, 0.2e1 * t54, pkin(3) ^ 2 + t37 * t48 + t33, t30, t23, t22, t28, t21, t37, -0.2e1 * t9 * t42 + 0.2e1 * t24, 0.2e1 * t59, t50 * t66, t10 ^ 2 + t9 ^ 2 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t40, 0, t16, 0, 0, 0, 0, 0, 0, -t56, t29, -t19, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t40, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t56, t29, -t19, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t57, -t42, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t57, -t42, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t44, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
