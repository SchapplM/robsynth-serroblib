% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:04
% EndTime: 2019-12-31 18:52:07
% DurationCPUTime: 0.68s
% Computational Cost: add. (543->75), mult. (1097->147), div. (0->0), fcn. (1229->6), ass. (0->62)
t41 = sin(pkin(8));
t42 = cos(pkin(8));
t44 = sin(qJ(3));
t64 = cos(qJ(3));
t28 = t41 * t64 + t44 * t42;
t70 = -0.2e1 * t28;
t59 = pkin(6) + qJ(2);
t29 = t59 * t42;
t54 = t59 * t41;
t12 = t29 * t44 + t64 * t54;
t69 = t12 ^ 2;
t26 = t44 * t41 - t42 * t64;
t23 = t26 ^ 2;
t35 = -pkin(2) * t42 - pkin(1);
t68 = 0.2e1 * t35;
t67 = 0.2e1 * t42;
t66 = t26 * pkin(4);
t45 = cos(qJ(4));
t65 = t45 * pkin(4);
t43 = sin(qJ(4));
t63 = t43 * t26;
t62 = t43 * t28;
t61 = t43 * t45;
t14 = t29 * t64 - t44 * t54;
t60 = t45 * t14;
t22 = t45 * t28;
t58 = -qJ(5) - pkin(7);
t37 = t41 ^ 2;
t38 = t42 ^ 2;
t57 = t37 + t38;
t39 = t43 ^ 2;
t40 = t45 ^ 2;
t32 = t39 + t40;
t56 = qJ(5) * t28;
t55 = t26 * t70;
t7 = pkin(3) * t26 - pkin(7) * t28 + t35;
t3 = -t14 * t43 + t45 * t7;
t53 = -pkin(3) * t28 - pkin(7) * t26;
t48 = -t45 * t56 + t3;
t1 = t48 + t66;
t2 = t60 + (t7 - t56) * t43;
t52 = t1 * t45 + t2 * t43;
t4 = t43 * t7 + t60;
t51 = t3 * t45 + t4 * t43;
t50 = -t3 * t43 + t4 * t45;
t30 = t58 * t43;
t31 = t58 * t45;
t49 = t30 * t45 - t31 * t43;
t36 = -pkin(3) - t65;
t33 = 0.2e1 * t61;
t24 = t28 ^ 2;
t21 = t45 * t26;
t20 = t40 * t24;
t17 = t39 * t24;
t16 = t43 * t22;
t15 = -0.2e1 * t24 * t61;
t11 = 0.2e1 * t26 * t22;
t10 = t43 * t55;
t9 = t32 * t28;
t8 = (-t39 + t40) * t28;
t5 = pkin(4) * t62 + t12;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, t41 * t67, 0, t38, 0, 0, pkin(1) * t67, -0.2e1 * pkin(1) * t41, 0.2e1 * t57 * qJ(2), qJ(2) ^ 2 * t57 + pkin(1) ^ 2, t24, t55, 0, t23, 0, 0, t26 * t68, t28 * t68, 0.2e1 * t12 * t28 - 0.2e1 * t14 * t26, t14 ^ 2 + t35 ^ 2 + t69, t20, t15, t11, t17, t10, t23, 0.2e1 * t12 * t62 + 0.2e1 * t26 * t3, 0.2e1 * t12 * t22 - 0.2e1 * t26 * t4, t51 * t70, t3 ^ 2 + t4 ^ 2 + t69, t20, t15, t11, t17, t10, t23, 0.2e1 * t1 * t26 + 0.2e1 * t5 * t62, -0.2e1 * t2 * t26 + 0.2e1 * t22 * t5, t52 * t70, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t41, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t26, t28, 0, t35, 0, 0, 0, 0, 0, 0, t21, -t63, -t9, t51, 0, 0, 0, 0, 0, 0, t21, -t63, -t9, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, 0, -t12, -t14, 0, 0, t16, t8, t63, -t16, t21, 0, -t12 * t45 + t43 * t53, t12 * t43 + t45 * t53, t50, -t12 * pkin(3) + pkin(7) * t50, t16, t8, t63, -t16, t21, 0, t26 * t30 + t36 * t62 - t45 * t5, t22 * t36 + t26 * t31 + t43 * t5, -t1 * t43 + t2 * t45 - t28 * t49, t1 * t30 - t2 * t31 + t36 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t39, t33, 0, t40, 0, 0, 0.2e1 * pkin(3) * t45, -0.2e1 * pkin(3) * t43, 0.2e1 * t32 * pkin(7), pkin(7) ^ 2 * t32 + pkin(3) ^ 2, t39, t33, 0, t40, 0, 0, -0.2e1 * t36 * t45, 0.2e1 * t36 * t43, -0.2e1 * t30 * t43 - 0.2e1 * t31 * t45, t30 ^ 2 + t31 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t62, t26, t3, -t4, 0, 0, 0, 0, t22, 0, -t62, t26, t48 + 0.2e1 * t66, -t2, -pkin(4) * t22, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t43, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t43, 0, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t45, 0, -t43 * pkin(7), -t45 * pkin(7), 0, 0, 0, 0, t43, 0, t45, 0, t30, t31, -t43 * pkin(4), t30 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t22, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t43, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
