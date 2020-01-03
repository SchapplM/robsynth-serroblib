% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR13_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t40 = cos(qJ(4));
t29 = -t40 * pkin(4) - pkin(3);
t72 = 0.2e1 * t29;
t38 = sin(qJ(3));
t71 = 0.2e1 * t38;
t70 = 0.2e1 * t40;
t69 = 2 * qJ(2);
t68 = -pkin(8) - pkin(7);
t36 = sin(qJ(5));
t67 = t36 * pkin(4);
t39 = cos(qJ(5));
t66 = t39 * pkin(4);
t41 = cos(qJ(3));
t21 = t38 * pkin(3) - t41 * pkin(7) + qJ(2);
t37 = sin(qJ(4));
t42 = -pkin(1) - pkin(6);
t56 = t40 * t42;
t50 = t38 * t56;
t5 = t50 + (-pkin(8) * t41 + t21) * t37;
t65 = t39 * t5;
t64 = t41 * pkin(3);
t19 = t36 * t40 + t39 * t37;
t63 = t19 * t38;
t62 = t37 * t38;
t61 = t37 * t40;
t60 = t37 * t41;
t59 = t37 * t42;
t58 = t38 * t42;
t57 = t40 * t38;
t28 = t40 * t41;
t55 = t41 * t19;
t54 = t41 * t38;
t53 = t41 * t42;
t31 = t37 ^ 2;
t33 = t40 ^ 2;
t52 = t31 + t33;
t32 = t38 ^ 2;
t34 = t41 ^ 2;
t26 = t32 + t34;
t51 = -0.2e1 * t54;
t49 = t37 * t28;
t15 = t40 * t21;
t4 = -pkin(8) * t28 + t15 + (pkin(4) - t59) * t38;
t1 = -t36 * t5 + t39 * t4;
t48 = t52 * t38;
t47 = -pkin(7) * t38 - t64;
t8 = -t37 * t58 + t15;
t9 = t37 * t21 + t50;
t46 = -t8 * t37 + t9 * t40;
t43 = qJ(2) ^ 2;
t35 = t42 ^ 2;
t30 = t34 * t35;
t23 = t68 * t40;
t22 = t68 * t37;
t20 = t26 * t42;
t17 = t36 * t37 - t39 * t40;
t16 = (pkin(4) * t37 - t42) * t41;
t14 = t39 * t28 - t36 * t60;
t13 = -t36 * t62 + t39 * t57;
t7 = t36 * t22 - t39 * t23;
t6 = t39 * t22 + t36 * t23;
t2 = t36 * t4 + t65;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t69, pkin(1) ^ 2 + t43, t34, t51, 0, t32, 0, 0, t38 * t69, t41 * t69, -0.2e1 * t20, t32 * t35 + t30 + t43, t33 * t34, -0.2e1 * t34 * t61, t54 * t70, t31 * t34, t37 * t51, t32, -0.2e1 * t34 * t59 + 0.2e1 * t8 * t38, -0.2e1 * t34 * t56 - 0.2e1 * t9 * t38, 0.2e1 * (-t37 * t9 - t40 * t8) * t41, t8 ^ 2 + t9 ^ 2 + t30, t14 ^ 2, -0.2e1 * t14 * t55, t14 * t71, t55 ^ 2, -t55 * t71, t32, 0.2e1 * t1 * t38 + 0.2e1 * t16 * t55, 0.2e1 * t16 * t14 - 0.2e1 * t2 * t38, -0.2e1 * t1 * t14 - 0.2e1 * t2 * t55, t1 ^ 2 + t16 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t26, t20, 0, 0, 0, 0, 0, 0, -t26 * t37, -t26 * t40, 0, t34 * t42 + t46 * t38, 0, 0, 0, 0, 0, 0, -t38 * t63 - t41 * t55, -t13 * t38 - t41 * t14, -t13 * t55 + t14 * t63, -t1 * t63 + t2 * t13 - t16 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t32 + t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 ^ 2 + t63 ^ 2 + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t38, 0, t53, -t58, 0, 0, t49, (-t31 + t33) * t41, t62, -t49, t57, 0, t47 * t37 + t40 * t53, -t37 * t53 + t47 * t40, t46, pkin(3) * t53 + t46 * pkin(7), t14 * t19, -t14 * t17 - t19 * t55, t63, t55 * t17, -t17 * t38, 0, t16 * t17 + t29 * t55 + t6 * t38, t29 * t14 + t16 * t19 - t7 * t38, -t1 * t19 - t6 * t14 - t2 * t17 - t55 * t7, t1 * t6 + t16 * t29 + t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t38, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t60, t48, pkin(7) * t48 + t64, 0, 0, 0, 0, 0, 0, -t41 * t17, -t55, -t13 * t17 + t19 * t63, t13 * t7 - t41 * t29 - t6 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t31, 0.2e1 * t61, 0, t33, 0, 0, pkin(3) * t70, -0.2e1 * pkin(3) * t37, 0.2e1 * t52 * pkin(7), t52 * pkin(7) ^ 2 + pkin(3) ^ 2, t19 ^ 2, -0.2e1 * t19 * t17, 0, t17 ^ 2, 0, 0, t17 * t72, t19 * t72, -0.2e1 * t7 * t17 - 0.2e1 * t6 * t19, t29 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t60, t38, t8, -t9, 0, 0, 0, 0, t14, 0, -t55, t38, t38 * t66 + t1, -t65 + (-t38 * pkin(4) - t4) * t36, (-t14 * t39 - t36 * t55) * pkin(4), (t1 * t39 + t2 * t36) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t57, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t13, 0, (t13 * t36 - t39 * t63) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t40, 0, -t37 * pkin(7), -t40 * pkin(7), 0, 0, 0, 0, t19, 0, -t17, 0, t6, -t7, (-t17 * t36 - t19 * t39) * pkin(4), (t36 * t7 + t39 * t6) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t66, -0.2e1 * t67, 0, (t36 ^ 2 + t39 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, -t55, t38, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, 0, t6, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t66, -t67, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
